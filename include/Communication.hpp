//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_COMM_HPP
#define SPEC_COMM_HPP

#include <numeric>
#include <mpi.h>
struct CommunicationDatatype {
    MPI_Datatype element_datatype;
    MPI_Datatype minimal_datatype;
    CommunicationDatatype(const MPI_Datatype &el, const MPI_Datatype &min) : element_datatype(el), minimal_datatype(min){}
    void free_datatypes() { MPI_Type_free(&element_datatype); MPI_Type_free(&minimal_datatype);}
};

template<class A>
inline std::vector<A> gather_elements_on(const std::vector<A> &local_el,
                                         const int dest_rank,
                                         const MPI_Datatype& sendtype,
                                         const MPI_Comm &comm) {
    std::vector<A> dest_el;
    int nlocal = local_el.size();
    int my_rank; MPI_Comm_rank(comm, &my_rank);
    int world_size; MPI_Comm_size(comm, &world_size);
    std::vector<int> counts(world_size, 0), displs(world_size, 0);
    MPI_Gather(&nlocal, 1, MPI_INT, &counts.front(), 1, MPI_INT, dest_rank, comm);
    int nb_elements = std::accumulate(counts.begin(), counts.end(), 0);
    for (int cpt = 0; cpt < world_size; ++cpt)
        displs[cpt] = cpt == 0 ? 0 : displs[cpt - 1] + counts[cpt - 1];
    if (my_rank == dest_rank) dest_el.resize(nb_elements);
    MPI_Gatherv(&local_el.front(), nlocal, sendtype,
                &dest_el.front(), &counts.front(), &displs.front(), sendtype, dest_rank, comm);

    return dest_el;
}

template<class A>
inline std::vector<A> gather_elements_on(const std::initializer_list<A> &&_local_el,
                                         const int dest_rank,
                                         const MPI_Datatype &sendtype,
                                         const MPI_Comm &comm) {
    std::vector<A> dest_el;
    std::vector<A> local_el(_local_el.begin(), _local_el.end());
    int nlocal = local_el.size();
    int my_rank; MPI_Comm_rank(comm, &my_rank);
    int world_size; MPI_Comm_size(comm, &world_size);
    std::vector<int> counts(world_size, 0), displs(world_size, 0);
    MPI_Gather(&nlocal, 1, MPI_INT, &counts.front(), 1, MPI_INT, dest_rank, comm);
    int nb_elements = std::accumulate(counts.begin(), counts.end(), 0);
    for (int cpt = 0; cpt < world_size; ++cpt) displs[cpt] = cpt == 0 ? 0 : displs[cpt - 1] + counts[cpt - 1];
    if (my_rank == dest_rank) dest_el.resize(nb_elements);
    MPI_Gatherv(&local_el.front(), nlocal, sendtype,
                &dest_el.front(), &counts.front(), &displs.front(), sendtype, dest_rank, comm);

    return dest_el;
}

template<class A>
inline void gather_elements_on(const std::vector<A> &local_el,
                               const int dest_rank,
                               std::vector<A>* _dest_el,
                               const MPI_Datatype &sendtype,
                               const MPI_Comm &comm) {
    std::vector<A>& dest_el = *(_dest_el), buffer;
    int nlocal = local_el.size();
    int my_rank; MPI_Comm_rank(comm, &my_rank);
    int world_size; MPI_Comm_size(comm, &world_size);
    std::vector<int> counts(world_size, 0), displs(world_size, 0);
    MPI_Gather(&nlocal, 1, MPI_INT, &counts.front(), 1, MPI_INT, dest_rank, comm);
    unsigned int nb_elements = std::accumulate(counts.begin(), counts.end(), 0);
    for (int cpt = 0; cpt < world_size; ++cpt) displs[cpt] = cpt == 0 ? 0 : displs[cpt - 1] + counts[cpt - 1];
    //if(_dest_el->size() == nb_elements)
    //if (my_rank == dest_rank) dest_el.resize(nb_elements);
    MPI_Gatherv(&local_el.front(), nlocal, sendtype,
                &_dest_el->front(), &counts.front(), &displs.front(), sendtype, dest_rank, comm);
    /*else {
        if (my_rank == dest_rank){
            buffer.resize(nb_elements);
            std::move(buffer.begin(), buffer.end(), std::back_inserter(dest_el));
        }
    }*/
}
#ifdef WITH_ZOLTAN
#include <zoltan.h>
template<class A>
std::tuple<const std::vector<A>, std::vector<std::vector<unsigned long>> >  zoltan_exchange_data(Zoltan_Struct *load_balancer,
                                          const std::vector<A> &data,
                                          int* nb_elements_recv,
                                          int* nb_elements_sent,
                                          MPI_Datatype datatype,
                                          MPI_Comm LB_COMM,
                                          double cell_size = 1.0) {
    *nb_elements_recv=0;
    *nb_elements_sent=0;

    int wsize;
    MPI_Comm_size(LB_COMM, &wsize);
    int caller_rank;
    MPI_Comm_rank(LB_COMM, &caller_rank);

    std::vector<A> buffer;
    std::vector<A> remote_data_gathered;
    long nb_reqs = 0;
    if (wsize == 1) return {remote_data_gathered, {}};

    std::vector<std::vector<A>> data_to_migrate(wsize);
    std::vector<std::vector<unsigned long>> id_to_migrate(wsize);
    size_t data_id = 0;
    std::vector<int> PEs(wsize, -1), parts(wsize, -1);
    int num_found, num_known = 0;
    std::vector<int> export_gids, export_lids, export_procs;

    // so much memory could be allocated here... potentially PE * n * DIM * 44 bytes => so linear in N
    // as DIM << PE <<<< n
    while (data_id < data.size()) {
        auto pos_in_double = data.at(data_id).template get_position_as_array<double>();
        int numparts;
        Zoltan_LB_Box_PP_Assign(load_balancer,
                             pos_in_double.at(0) - cell_size,
                             pos_in_double.at(1) - cell_size,
                             pos_in_double.size() == 3 ? pos_in_double.at(2) - cell_size : 0.0,
                             pos_in_double.at(0) + cell_size,
                             pos_in_double.at(1) + cell_size,
                             pos_in_double.size() == 3 ? pos_in_double.at(2) + cell_size : 0.0,
                             &PEs.front(), &num_found, parts.data(), &numparts);

        for (int PE_idx = 0; PE_idx < num_found; PE_idx++) {
            int PE = PEs[PE_idx];
            if (PE >= 0) {
                if (PE != caller_rank) {
                    export_gids.push_back(data.at(data_id).gid);
                    export_lids.push_back(data_id);
                    export_procs.push_back(PE);
                    //get the value and copy it into the "to migrate" vector
                    data_to_migrate.at(PE).push_back(data.at(data_id));
                    id_to_migrate.at(PE).push_back(data_id);
                    num_known++;
                }
            }
        }
        data_id++; //if the element must stay with me then check the next one
    }
// Compute who has to send me something via Zoltan.
///    int ierr = Zoltan_Invert_Lists(load_balancer, num_known, known_gids, known_lids, &export_procs[0], &export_procs[0],
///                                   &num_found, &found_gids, &found_lids, &found_procs, &found_parts);

    std::vector<int> num_import_from_procs(wsize);
    std::vector<int> import_from_procs;

    {
        nb_reqs = std::count_if(data_to_migrate.cbegin(), data_to_migrate.cend(), [](auto buf){return !buf.empty();});
        std::vector<MPI_Request> send_reqs(nb_reqs),
                rcv_reqs(wsize);
        std::vector<MPI_Status> statuses(wsize);

        int nb_neighbor = 0;
        for(int comm_pe = 0; comm_pe < wsize; ++comm_pe){
            MPI_Irecv(&num_import_from_procs[comm_pe], 1, MPI_INT, comm_pe, 666, LB_COMM, &rcv_reqs[comm_pe]);
            int send_size = data_to_migrate.at(comm_pe).size();
            if (send_size) {
                MPI_Isend(&send_size, 1, MPI_INT, comm_pe, 666, LB_COMM, &send_reqs[nb_neighbor]);
                nb_neighbor++;
            }
        }
        MPI_Waitall(send_reqs.size(), &send_reqs.front(), MPI_STATUSES_IGNORE);
        MPI_Barrier(LB_COMM);
        for(int comm_pe = 0; comm_pe < wsize; ++comm_pe) {
            int flag; MPI_Status status;
            MPI_Test(&rcv_reqs[comm_pe], &flag, &status);
            if(!flag) MPI_Cancel(&rcv_reqs[comm_pe]);
        }

        import_from_procs.clear();
        for(int i = 0; i < wsize; ++i){
            if(num_import_from_procs[i] > 0) import_from_procs.push_back(i);
        }
    }

    int cpt = 0;

// Send the data to neighbors
    std::vector<MPI_Request> reqs(nb_reqs);
    *nb_elements_sent = 0;
    for (size_t PE = 0; PE < (unsigned int) wsize; PE++) {
        int send_size = data_to_migrate.at(PE).size();
        if (send_size) {
            *nb_elements_sent += send_size;
            MPI_Isend(&data_to_migrate.at(PE).front(), send_size, datatype, PE, 400, LB_COMM,
                      &reqs[cpt]);
            cpt++;
        }
    }
// Import the data from neighbors
    *nb_elements_recv = 0;
    for (int proc_id : import_from_procs) {
        size_t size = num_import_from_procs[proc_id];
        nb_elements_recv += size;
        buffer.resize(size);
        MPI_Recv(&buffer.front(), size, datatype, proc_id, 400, LB_COMM, MPI_STATUS_IGNORE);
        std::move(buffer.begin(), buffer.end(), std::back_inserter(remote_data_gathered));
    }

    MPI_Waitall(reqs.size(), &reqs.front(), MPI_STATUSES_IGNORE);
    return std::make_tuple(remote_data_gathered, id_to_migrate);
}

template<class A>
void zoltan_migrate_particles(Zoltan_Struct *load_balancer,
                              std::vector<A> *_data,
                              const MPI_Datatype datatype,
                              const MPI_Comm LB_COMM) {
    std::vector<A>& data = *(_data);
    int wsize;
    MPI_Comm_size(LB_COMM, &wsize);
    int caller_rank;
    MPI_Comm_rank(LB_COMM, &caller_rank);

    if (wsize == 1) return;

    std::vector<std::vector<A>> data_to_migrate(wsize);

    size_t data_id = 0;
    int PE;
    int num_known = 0;
    std::vector<int> export_gids, export_lids, export_procs;
    while (data_id < data.size()) {
        auto pos_in_double = data.at(data_id).template get_position_as_array<double>();
        Zoltan_LB_Point_Assign(load_balancer, &pos_in_double.front(), &PE);
        if (PE != caller_rank) {
            export_gids.push_back(data.at(data_id).gid);
            export_lids.push_back(data.at(data_id).lid);
            export_procs.push_back(PE);
            //if the current element has to be moved, then swap with the last and pop it out (dont need to move the pointer also)
            //swap iterator values in constant time
            std::iter_swap(data.begin()+ data_id, data.end() - 1);
// get the value and push it in the "to migrate" vector
            data_to_migrate.at(PE).push_back(*(data.end()- 1));
// pop the head of the list in constant time
            data.pop_back();
            num_known++;
        } else data_id++; //if the element must stay with me then check the next one
    }

    ZOLTAN_ID_PTR known_gids = (ZOLTAN_ID_PTR) &export_gids.front();
    ZOLTAN_ID_PTR known_lids = (ZOLTAN_ID_PTR) &export_lids.front();
    ZOLTAN_ID_PTR found_gids, found_lids;

    int *found_procs, *found_parts, num_found;

    int ierr = Zoltan_Invert_Lists(load_balancer, num_known, known_gids, known_lids, &export_procs[0], &export_procs[0],
                                   &num_found, &found_gids, &found_lids, &found_procs, &found_parts);

    std::vector<int> num_import_from_procs(wsize);
    std::vector<int> import_from_procs;

    for (size_t i = 0;i < num_found;++i) {
        num_import_from_procs[found_procs[i]]++;
        if (std::find(import_from_procs.begin(), import_from_procs.end(), found_procs[i]) == import_from_procs.end())
            import_from_procs.push_back(found_procs[i]);
    }

    /* Let's Migrate ma boi ! */

    if (num_found > 0)
        Zoltan_LB_Free_Part(&found_gids, &found_lids, &found_procs, &found_parts);

    int nb_reqs = 0;
    for (auto buf: data_to_migrate) {
        if (!buf.empty())
            nb_reqs++;
    }

    int cpt = 0;
    std::vector<MPI_Request> reqs(nb_reqs);
    for (size_t PE = 0;PE < wsize;PE++) {
        int send_size = data_to_migrate.at(PE).size();
        if (send_size) {
            MPI_Isend(&data_to_migrate.at(PE).front(), send_size, datatype, PE, 300, LB_COMM, &reqs[cpt]);
            cpt++;
        }
    }
    std::vector<A> buffer;
    for (int proc_id : import_from_procs) {
        size_t size = num_import_from_procs[proc_id];
        buffer.resize(size);
        MPI_Recv(&buffer.front(), size, datatype, proc_id, 300, LB_COMM, MPI_STATUS_IGNORE);
        std::move(buffer.begin(), buffer.end(), std::back_inserter(data));
    }

    const auto nb_data = data.size();

    for (int i = 0;i < nb_data;++i)
        data[i].lid = i;

    MPI_Waitall(reqs.size(), &reqs.front(), MPI_STATUSES_IGNORE);

}

inline bool point_belongs_to_me(Zoltan_Struct* load_balancer, std::array<double, 2> pos, int my_rank){
    int PE;
    Zoltan_LB_Point_Assign(load_balancer, &pos.front(), &PE);
    return PE == my_rank;
}
#endif

template<class A>
std::tuple<long, std::vector<int>, std::vector<int>> compute_invert_list(const std::vector<std::vector<A>>& data_to_migrate, const int worldsize, MPI_Comm& world) {

    long nb_reqs;
    std::vector<int> num_import_from_procs(worldsize), import_from_procs;

    nb_reqs = std::count_if(data_to_migrate.cbegin(), data_to_migrate.cend(), [](auto buf){return !buf.empty();});
    std::vector<MPI_Request> send_reqs(nb_reqs), rcv_reqs(worldsize);
    std::vector<MPI_Status> statuses(worldsize);

    int nb_neighbor = 0;
    for(int comm_pe = 0; comm_pe < worldsize; ++comm_pe){
        MPI_Irecv(&num_import_from_procs[comm_pe], 1, MPI_INT, comm_pe, 666, world, &rcv_reqs[comm_pe]);
        int send_size = data_to_migrate.at(comm_pe).size();
        if (send_size) {
            MPI_Isend(&send_size, 1, MPI_INT, comm_pe, 666, world, &send_reqs[nb_neighbor]);
            nb_neighbor++;
        }
    }
    MPI_Waitall(send_reqs.size(), &send_reqs.front(), MPI_STATUSES_IGNORE);
    MPI_Barrier(world);
    for(int comm_pe = 0; comm_pe < worldsize; ++comm_pe) {
        int flag; MPI_Status status;
        MPI_Test(&rcv_reqs[comm_pe], &flag, &status);
        if(!flag) MPI_Cancel(&rcv_reqs[comm_pe]);
    }

    import_from_procs.clear();
    for(int i = 0; i < worldsize; ++i){
        if(num_import_from_procs[i] > 0) import_from_procs.push_back(i);
    }

    return std::make_tuple(nb_reqs, num_import_from_procs, import_from_procs);
}


#endif //SPEC_COMM_HPP
