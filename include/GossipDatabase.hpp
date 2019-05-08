//
// Created by xetql on 2/23/19.
//

#ifndef SPEC_CPULOADDATABASE_HPP
#define SPEC_CPULOADDATABASE_HPP

#include <cassert>
#include <tuple>
#include <vector>
#include <random>
#include <memory>
#include <mpi.h>
#include <algorithm>
#include <ostream>
#include <queue>
#include "Utils.hpp"
#include <type_traits>
#include <functional>


template<class StoredDataType>
class GossipDatabase {
    typedef int Index;
    typedef int Age;
    struct DatabaseEntry {
        using EntryUpdateStrategy = std::function<DatabaseEntry (DatabaseEntry, DatabaseEntry)>;
        Index idx;
        Age age;
        StoredDataType load;

        DatabaseEntry() :
                idx(-1),
                age(std::numeric_limits<int>::max()),
                load(-1)
        {};

        DatabaseEntry(Index _idx, Age _age, StoredDataType _load) :
                idx(_idx),
                age(_age),
                load(_load)
        {};

        DatabaseEntry(Index _idx, StoredDataType _load) :
                idx(_idx),
                age(0),
                load(_load)
        {};

        friend std::ostream &operator<<(std::ostream &os, const DatabaseEntry &entry) {
            os << "idx: " << entry.idx << " age: " << entry.age << " load: " << entry.load;
            return os;
        }
    };
    using EntryPropagationPredicate = std::function<bool(DatabaseEntry)>;
    EntryPropagationPredicate default_predicate = [](auto e) { return e.idx >= 0; };

    using EntryUpdateStrategy = std::function<DatabaseEntry (DatabaseEntry, DatabaseEntry)>;
    EntryUpdateStrategy update_strategy = GossipDatabase::default_update_strategy; //[&](DatabaseEntry old, DatabaseEntry mine) { };

    using value_type = typename std::enable_if<
            std::is_same<StoredDataType , double>::value ||
            std::is_same<StoredDataType , float>::value  ||
            std::is_same<StoredDataType , long>::value   ||
            std::is_same<StoredDataType , int>::value   >::type;

    MPI_Datatype entry_datatype;
    MPI_Status current_recv_status;
    const int database_size, number_of_message = 2, SEND_TAG;
    MPI_Comm world;
    std::vector<DatabaseEntry> pe_load_data;

    int my_rank, worldsize;

    mutable std::vector<MPI_Request> current_send_reqs;
    mutable int  step_counter = 0;
    mutable std::vector<std::vector<DatabaseEntry>> snd_entries;

    static DatabaseEntry default_update_strategy(DatabaseEntry old, DatabaseEntry mine) {
        return old.idx == mine.idx ? mine : is_initialized(old) ? DatabaseEntry(old.idx, old.age+1, old.load) : old;
    }

public:

    GossipDatabase(int database_size, int number_of_message, int send_tag, MPI_Comm _world, EntryUpdateStrategy strategy) :
            update_strategy(strategy),
            database_size(database_size),
            number_of_message(number_of_message),
            SEND_TAG(send_tag),
            world(_world) {
        MPI_Comm_rank(world, &my_rank);
        MPI_Comm_size(world, &worldsize);
        for(int i = 0; i < database_size; ++i)
            pe_load_data.emplace_back(i, std::numeric_limits<int>::max(), 0.0);
        register_datatype();
        srand(my_rank + time(NULL));
        snd_entries.resize(number_of_message);
        current_send_reqs.resize(number_of_message);
    }

    GossipDatabase(int database_size, int number_of_message, int send_tag, MPI_Comm _world) :
        database_size(database_size),
        number_of_message(number_of_message),
        SEND_TAG(send_tag),
        world(_world) {
        MPI_Comm_rank(world, &my_rank);
        MPI_Comm_size(world, &worldsize);
        for(int i = 0; i < database_size; ++i)
            pe_load_data.emplace_back(i, std::numeric_limits<int>::max(), 0.0);
        register_datatype();
        srand(my_rank + time(NULL));
        snd_entries.resize(number_of_message);
        current_send_reqs.resize(number_of_message);
    }

    void execute(Index idx, StoredDataType data) {
        execute(idx, data, update_strategy, default_predicate);
    }

    void execute(Index idx, StoredDataType data, EntryUpdateStrategy strategy, EntryPropagationPredicate pred) {
        if(worldsize > 2){
            if(step_counter > 0){
                finish_gossip_step();
            }
            gossip_update(idx, data, strategy);
            gossip_propagate(pred);
        }
    }

    StoredDataType get(Index idx) const {
        return pe_load_data[idx].load;
    }

    inline void reset(){
        pe_load_data.clear();
        pe_load_data.resize(database_size);
        step_counter = 0;
    }



    StoredDataType skewness() const {
        std::vector<StoredDataType> &&loads = get_all_meaningfull_data();
        return stats::skewness<StoredDataType>(loads.begin(), loads.end());
    }

    StoredDataType mean() const {
        std::vector<StoredDataType> &&loads = get_all_meaningfull_data();
        return stats::mean<StoredDataType>(loads.begin(), loads.end());
    }

    int max_age() const {
        int age = -1;
        for(const auto& entry : pe_load_data) age = std::max(age, entry.age);
        return age == -1 ? std::numeric_limits<int>::max() : age;
    }

    bool has_converged(Age threshold) const {
        return std::all_of(pe_load_data.begin(), pe_load_data.end(), [&threshold](auto entry){return entry.age < threshold;});
    }

    StoredDataType sum() const {
        return this->mean() * worldsize;
    }

    StoredDataType variance() const {
        const auto mu = mean();
        auto data = get_all_meaningfull_data();
        auto N = data.size();
        StoredDataType var = 0.0;

        for(auto& x : data) var += std::pow(x-mu, 2.0);
        var /= (StoredDataType) N;
        //var -= std::pow(mu, 2.0);

        return var;
    }

    double zscore(Index idx) const {
        const auto load = get(idx);
        const auto mu = mean();
        const auto stddev = std::sqrt(variance());
        return (load - mu) / stddev;
    }

    std::vector<StoredDataType> get_all_data() const {
        std::vector<StoredDataType> loads(worldsize, -1);
        for (auto it = pe_load_data.begin(); it != pe_load_data.end(); it++) {
            //auto i = std::distance(pe_load_data.begin(), it);
            auto& entry = *it;
            if(entry.idx >= 0) loads[entry.idx] = ((*it).load);
        }
        return loads;
    }

private:
    void free_datatypes() {
        MPI_Type_free(&entry_datatype);
    }

    void set(Index idx, Age age, StoredDataType my_load) {
        pe_load_data[idx] = {idx, age, my_load};
    }

    /**
     * Randomly select a processing elements and "contaminate" him with my information
     */
    void gossip_propagate() {
        gossip_propagate([](auto e) { return e.idx >= 0; });
    }

    void gossip_propagate(EntryPropagationPredicate pred) const {
        std::vector<int> destinations;
        int destination;
        for(int i = 0; i < number_of_message; ++i) {
            do {
                destination = rand() % worldsize;
            } while(std::find(destinations.begin(), destinations.end(), destination) != destinations.end() || destination == my_rank);
            destinations.push_back(destination);
            snd_entries[i].clear();
            // propagate entries that match a predicate
            std::copy_if(pe_load_data.begin(), pe_load_data.end(), std::back_inserter(snd_entries[i]), pred);
            MPI_Isend(snd_entries[i].data(), snd_entries[i].size(), entry_datatype, destination, SEND_TAG, world,
                      &current_send_reqs[i]);
        }
        step_counter++;
    }

    /**
     * Update information ages and my load
     */
    void gossip_update(Index idx, StoredDataType my_load) {
        gossip_update(idx, 0, my_load, update_strategy);
    }

    void gossip_update(Index idx, StoredDataType my_load, EntryUpdateStrategy strategy) {
        assert(idx < database_size);
        assert(idx >= 0);
        DatabaseEntry mine = {idx, 0, my_load};
        for (DatabaseEntry &entry : pe_load_data) entry = strategy(entry, mine);
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

    void register_datatype() {
        MPI_Datatype oldtype_element[2];
        MPI_Aint offset[2], lb, int_offset;
        const int number_of_int_elements = 2;
        const int number_of_data_elements = 1;
        int blockcount_element[2];
        blockcount_element[0] = number_of_int_elements; // gid, lid, exit, waiting_time
        blockcount_element[1] = number_of_data_elements; // position <x,y>
        oldtype_element[0] = MPI_INT;

        if(std::is_same<StoredDataType, int>::value) oldtype_element[1]         = MPI_INT;
        else if(std::is_same<StoredDataType, long>::value) oldtype_element[1]   = MPI_LONG;
        else if(std::is_same<StoredDataType, float>::value) oldtype_element[1]  = MPI_FLOAT;
        else if(std::is_same<StoredDataType, double>::value) oldtype_element[1] = MPI_DOUBLE;

        MPI_Type_get_extent(MPI_INT, &lb, &int_offset);
        offset[0] = static_cast<MPI_Aint>(0);
        offset[1] = number_of_int_elements * int_offset;
        MPI_Type_create_struct(2, blockcount_element, offset, oldtype_element, &entry_datatype);
        MPI_Type_commit(&entry_datatype);
    }

    std::vector<StoredDataType> get_all_meaningfull_data() const {
        std::vector<StoredDataType> loads;
        for (auto it = pe_load_data.begin(); it != pe_load_data.end(); it++) {
            //auto i = std::distance(pe_load_data.begin(), it);
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

    static bool is_initialized(const DatabaseEntry& entry) {
        return entry.age < std::numeric_limits<int>::max();
    }

};


#endif //SPEC_CPULOADDATABASE_HPP
