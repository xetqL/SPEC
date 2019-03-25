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
#include <queue>
class CPULoadDatabase {
    typedef int Index;
    typedef int Age;
    typedef double PELoad;

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
    mutable std::vector<MPI_Request> current_send_reqs;
    MPI_Status current_recv_status;

    const int database_size, number_of_message = 2, SEND_TAG;
    int current_recv_flag;
    int my_rank, worldsize;
    std::vector<DatabaseEntry> pe_load_data;
    mutable std::vector<std::vector<DatabaseEntry>> snd_entries;
public:

    CPULoadDatabase(int database_size, int number_of_message, int send_tag, MPI_Comm _world) :
        database_size(database_size),
        number_of_message(number_of_message),
        SEND_TAG(send_tag),
        world(_world)
    {
        MPI_Comm_rank(world, &my_rank);
        MPI_Comm_size(world, &worldsize);
        for(int i = 0; i < database_size; ++i){
            pe_load_data.emplace_back(i, std::numeric_limits<int>::max(), 0.0);
        }
        register_datatype();
        srand(my_rank + time(NULL));
        snd_entries.resize(number_of_message);
        current_send_reqs.resize(number_of_message);
    }

    void free_datatypes() {
        MPI_Type_free(&entry_datatype);
    }

    PELoad get(Index idx) const {
        return pe_load_data[idx].load;
    }

    void set(Index idx, Age age, PELoad my_load) {
        pe_load_data[idx] = {idx, age, my_load};
    }

    inline void reset(){
        pe_load_data.clear();
        pe_load_data.resize(database_size);
    }

    /**
     * Update information ages and my load
     */
    void gossip_update(Index idx, PELoad my_load) {
        gossip_update(idx, 0, my_load, [&](auto old, auto mine){ return old.idx == mine.idx ? mine : is_initialized(old) ? DatabaseEntry(old.idx, old.age+1, old.load) : old;});
    }

    template<class Strategy>
    void gossip_update(Index idx, Age age, PELoad my_load, Strategy strategy) {
        assert(idx < database_size);
        assert(idx >= 0);
        DatabaseEntry mine = {idx, age, my_load};
        for (DatabaseEntry &entry : pe_load_data) {
            entry = strategy(entry, mine);
        }
    }

    /**
     * Randomly select a processing elements and "contaminate" him with my information
     */
    void gossip_propagate() {
        gossip_propagate([](auto e) { return e.idx >= 0; });
    }

    template<class Predicate>
    void gossip_propagate(Predicate pred) const {
        std::vector<int> destinations;
        int destination;
        for(int i = 0; i < number_of_message; ++i) {
            do {
                destination = rand() % worldsize;
            } while(std::find(destinations.begin(), destinations.end(), destination) != destinations.end() || destination == my_rank);
            destinations.push_back(destination);
            snd_entries[i].clear();
            std::copy_if(pe_load_data.begin(), pe_load_data.end(), std::back_inserter(snd_entries[i]), pred);
            MPI_Isend(snd_entries[i].data(), snd_entries[i].size(), entry_datatype, destination, SEND_TAG, world,
                      &current_send_reqs[i]);
        }
    }

    void finish_gossip_step() {
        int more_messages = 1;
        while (more_messages) {
            MPI_Iprobe(MPI_ANY_SOURCE, SEND_TAG, world, &more_messages, &current_recv_status);
            if (more_messages) {
                int cnt;
                MPI_Get_count(&current_recv_status, entry_datatype, &cnt);
                std::vector<DatabaseEntry> rcv_entries(cnt);
                MPI_Recv(&rcv_entries.front(), cnt, entry_datatype, current_recv_status.MPI_SOURCE,
                         SEND_TAG, world, MPI_STATUS_IGNORE);
                merge_into_database(std::move(rcv_entries));
            }
        }
        MPI_Waitall(number_of_message, current_send_reqs.data(), MPI_STATUSES_IGNORE);
    }

    PELoad skewness() const {
        std::vector<PELoad> &&loads = get_all_meaningfull_data();
        return stats::skewness<PELoad>(loads.begin(), loads.end());
    }

    PELoad mean() const {
        std::vector<PELoad> &&loads = get_all_meaningfull_data();
        return stats::mean<PELoad>(loads.begin(), loads.end());
    }

    int max_age() const {
        int age = -1;
        for(const auto& entry : pe_load_data) age = std::max(age, entry.age);
        return age == -1 ? std::numeric_limits<int>::max() : age;
    }

    bool has_converged(Age threshold) const {
        return std::all_of(pe_load_data.begin(), pe_load_data.end(), [&threshold](auto entry){return entry.age < threshold;});
    }

    PELoad sum() const {
        return this->mean() * worldsize;
    }

    PELoad variance() const {
        const auto mu = mean();
        auto data = get_all_meaningfull_data();
        auto N = data.size();
        PELoad var = 0.0;

        for(auto& x : data) var += std::pow(x-mu, 2.0);
        var /= (PELoad) N;
        //var -= std::pow(mu, 2.0);

        return var;
    }

    double zscore(Index idx) const {
        const auto load = get(idx);
        const auto mu = mean();
        const auto stddev = std::sqrt(variance());
        return (load - mu) / stddev;
    }

    std::vector<PELoad> get_all_data() const {
        std::vector<PELoad> loads(worldsize, -1);
        for (auto it = pe_load_data.begin(); it != pe_load_data.end(); it++) {
            auto i = std::distance(pe_load_data.begin(), it);
            auto& entry = *it;
            if(entry.idx >= 0) loads[entry.idx] = ((*it).load);
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

    std::vector<PELoad> get_all_meaningfull_data() const {
        std::vector<PELoad> loads;
        for (auto it = pe_load_data.begin(); it != pe_load_data.end(); it++) {
            auto i = std::distance(pe_load_data.begin(), it);
            auto& entry = *it;
            if(entry.idx >= 0)
                loads.push_back((*it).load);
        }
        return loads;
    }

    template<class Strategy>
    void merge_into_database(std::vector<DatabaseEntry> &&data, Strategy&& strategy) {
        for (auto new_entry : data) {
            auto &entry = pe_load_data[new_entry.idx];
            entry = strategy(entry, new_entry);
        }
    }

    /**
     * Merge by keeping the most recent data
     */
    void merge_into_database(std::vector<DatabaseEntry> &&data) {
        merge_into_database(std::move(data), [](auto old, auto recv){ return old.age < recv.age ? old : recv;});
    }

    inline bool is_initialized(const DatabaseEntry& entry) {
        return entry.age < std::numeric_limits<int>::max();
    }

};


#endif //SPEC_CPULOADDATABASE_HPP
