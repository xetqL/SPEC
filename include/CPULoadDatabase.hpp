//
// Created by xetql on 2/23/19.
//

#ifndef SPEC_CPULOADDATABASE_HPP
#define SPEC_CPULOADDATABASE_HPP

#include <tuple>
#include <vector>
#include <random>
#include <memory>
#include <mpi.h>
#include <algorithm>
#include <ostream>

class CPULoadDatabase {
    typedef int Index;
    typedef int Age;
    typedef double PELoad;
    typedef double PESlope;

    const int SEND_TAG;

    struct DatabaseEntry {
        Index idx;
        Age age;
        PELoad load;

        DatabaseEntry() :
                idx(-1),
                age(std::numeric_limits<int>::max()),
                load(-1)
        {};

        DatabaseEntry(Index _idx, Age _age, PELoad _load) :
                idx(_idx),
                age(_age),
                load(_load)
        {};

        DatabaseEntry(Index _idx, PELoad _load) :
                idx(_idx),
                age(0),
                load(_load)
        {};

        friend std::ostream &operator<<(std::ostream &os, const DatabaseEntry &entry) {
            os << "idx: " << entry.idx << " age: " << entry.age << " load: " << entry.load;
            return os;
        }
    };

    MPI_Comm world;
    MPI_Datatype entry_datatype;
    MPI_Request current_send_reqs[2];
    MPI_Status current_recv_status;
    int current_recv_flag;
    int my_rank, worldsize;
    std::vector<DatabaseEntry> pe_load_data;
public:

    CPULoadDatabase(int send_tag, MPI_Comm _world) : SEND_TAG(send_tag), world(_world) {
        MPI_Comm_rank(world, &my_rank);
        MPI_Comm_size(world, &worldsize);
        pe_load_data.resize(worldsize);
        pe_load_data[my_rank].idx = my_rank;
        register_datatype();
        srand(my_rank + time(NULL));
    }
    void free_datatypes(){
        MPI_Type_free(&entry_datatype);
    }
    PELoad get(Index idx){
        return pe_load_data[idx].load;
    }

    inline void reset(){
        pe_load_data.clear(); pe_load_data.resize(worldsize);
    }

    /**
     * Update information ages and my load
     */
    void gossip_update(PELoad my_load) {
        for (DatabaseEntry &entry : pe_load_data) {
            if (entry.idx == my_rank) {
                entry.age = 0;
                entry.load = my_load;
                //entry.slope = my_slope;
            } else if (entry.idx >= 0) {
                entry.age++;
            }
        }
    }

    /**
     * Randomly select a processing elements and "contaminate" him with my information
     */
    void gossip_propagate() {
        int destination1, destination2;

        do {
            destination1 = rand() % worldsize;
        } while(destination1 == my_rank);

        do {
            destination2 = rand() % worldsize;
        } while(destination2 == destination1 || destination2 == my_rank);

        std::vector<DatabaseEntry> snd_entryA, snd_entryB;
        std::copy_if(pe_load_data.begin(), pe_load_data.end(), std::back_inserter(snd_entryA),
                [](auto e) { return e.idx >= 0; });
        MPI_Isend(&snd_entryA.front(), snd_entryA.size(), entry_datatype, destination1, SEND_TAG, world,
                  &current_send_reqs[0]);
        std::copy_if(pe_load_data.begin(), pe_load_data.end(), std::back_inserter(snd_entryB),
                     [](auto e) { return e.idx >= 0; });
        MPI_Isend(&snd_entryB.front(), snd_entryB.size(), entry_datatype, destination2, SEND_TAG, world,
                  &current_send_reqs[1]);
    }

    void finish_gossip_step() {

        MPI_Iprobe(MPI_ANY_SOURCE, SEND_TAG, world, &current_recv_flag, &current_recv_status);
        if (current_recv_flag) {
            int cnt;
            MPI_Get_count(&current_recv_status, entry_datatype, &cnt);
            std::vector<DatabaseEntry> rcv_entries(cnt);
            MPI_Recv(&rcv_entries.front(), cnt, entry_datatype, current_recv_status.MPI_SOURCE,
                     SEND_TAG, world, MPI_STATUS_IGNORE);
            merge_into_database(std::move(rcv_entries));
        }

        MPI_Iprobe(MPI_ANY_SOURCE, SEND_TAG, world, &current_recv_flag, &current_recv_status);
        if (current_recv_flag) {
            int cnt;
            MPI_Get_count(&current_recv_status, entry_datatype, &cnt);
            std::vector<DatabaseEntry> rcv_entries(cnt);
            MPI_Recv(&rcv_entries.front(), cnt, entry_datatype, current_recv_status.MPI_SOURCE,
                     SEND_TAG, world, MPI_STATUS_IGNORE);
            merge_into_database(std::move(rcv_entries));
        }

        MPI_Waitall(2, current_send_reqs, MPI_STATUSES_IGNORE);
        //if(!my_rank) std::cout << "Gossip done." << std::endl;

    }

    PELoad skewness() {
        std::vector<PELoad> &&loads = get_all_meaningfull_data();
        return stats::skewness<PELoad>(loads.begin(), loads.end());
    }

    PELoad mean() {
        std::vector<PELoad> &&loads = get_all_meaningfull_data();
        return stats::mean<PELoad>(loads.begin(), loads.end());
    }

    int max_age() {
        int age = -1;
        for(const auto& entry : pe_load_data){
            age = std::max(age, entry.age);
        }
        return age == -1 ? std::numeric_limits<int>::max() : age;
    }

    bool has_converged(Age threshold) {
        return std::all_of(pe_load_data.begin(), pe_load_data.end(), [&threshold](auto entry){return entry.age < threshold;});
    }

    PELoad sum() {
        return this->mean() * worldsize;
    }

    PELoad variance() {
        const auto mu = mean();
        auto data = get_all_meaningfull_data();
        auto N = data.size();
        PELoad var = 0.0;

        for(auto& x : data) var += std::pow(x-mu, 2.0);
        var /= (PELoad) N;
        //var -= std::pow(mu, 2.0);

        return var;
    }

    double zscore(Index idx) {
        const auto load = get(idx);
        const auto mu = mean();
        const auto stddev = std::sqrt(variance());
        return (load - mu) / stddev;
    }

    std::vector<PELoad> get_all_data() {
        std::vector<PELoad> loads(worldsize, -1);
        for (auto it = pe_load_data.begin(); it != pe_load_data.end(); it++) {
            auto i = std::distance(pe_load_data.begin(), it);
            auto& entry = *it;
            if(entry.idx >= 0)
                loads[entry.idx] = ((*it).load);
        }
        return loads;
    }

private:
    void register_datatype() {
        MPI_Datatype element_datatype,
                oldtype_element[2];
        MPI_Aint offset[2], intex, int_offset;
        const int number_of_int_elements = 2;
        const int number_of_double_elements = 1;
        int blockcount_element[2];
        blockcount_element[0] = number_of_int_elements; // gid, lid, exit, waiting_time
        blockcount_element[1] = number_of_double_elements; // position <x,y>
        oldtype_element[0] = MPI_INT;
        oldtype_element[1] = MPI_DOUBLE;
        MPI_Type_extent(MPI_INT, &int_offset);
        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = number_of_int_elements * int_offset;
        MPI_Type_struct(2, blockcount_element, offset, oldtype_element, &entry_datatype);
        MPI_Type_commit(&entry_datatype);
    }

    std::vector<PELoad> get_all_meaningfull_data() {
        std::vector<PELoad> loads;
        for (auto it = pe_load_data.begin(); it != pe_load_data.end(); it++) {
            auto i = std::distance(pe_load_data.begin(), it);
            auto& entry = *it;
            if(entry.idx >= 0)
                loads.push_back((*it).load);
        }
        return loads;
    }

    /**
     * Merge by keeping the most recent data
     */
    void merge_into_database(std::vector<DatabaseEntry> &&data) {
        for (auto new_entry : data) {
            auto &entry = pe_load_data[new_entry.idx];
            if (entry.age > new_entry.age) {
                entry = new_entry;
            }
        }
    }
};


#endif //SPEC_CPULOADDATABASE_HPP
