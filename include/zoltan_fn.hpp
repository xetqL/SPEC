//
// Created by xetql on 2/11/19.
//

#ifndef SPEC_ZOLTAN_FN_HPP
#define SPEC_ZOLTAN_FN_HPP
#include <cassert>
#include <random>
#include <string>
#include <vector>
#include <zoltan.h>
#include "communication.hpp"
#include "cell.hpp"
#include "utils.hpp"

#define ENABLE_AUTOMATIC_MIGRATION true
#define DISABLE_AUTOMATIC_MIGRATION FALSE


int get_number_of_objects(void *data, int *ierr) {
    std::vector<Cell> *mesh = (std::vector<Cell> *)data;
    *ierr = ZOLTAN_OK;
    return mesh->size();
}

void get_object_list(void *data, int sizeGID, int sizeLID,
                     ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                     int wgt_dim, float *obj_wgts, int *ierr) {
    size_t i;
    std::vector<Cell> *mesh = (std::vector<Cell> *)data;
    *ierr = ZOLTAN_OK;
    /* In this example, return the IDs of our objects, but no weights.
     * Zoltan will assume equally weighted objects.
     */
    for (i=0; i < mesh->size(); i++){
        globalID[i] = mesh->at(i).gid;
        localID[i] = i;
    }
}

int get_num_geometry(void *data, int *ierr) {
    *ierr = ZOLTAN_OK;
    return 2;
}

void get_geometry_list(void *data, int sizeGID, int sizeLID,
                       int num_obj,
                       ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int num_dim, double *geom_vec, int *ierr) {
    int i;

    std::vector<Cell> *mesh = (std::vector<Cell> *)data;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 2)){
        *ierr = ZOLTAN_FATAL;
        return;
    }

    *ierr = ZOLTAN_OK;

    for (i=0;  i < num_obj; i++){
        auto gid = mesh->at(i).gid;
        int x,y; std::tie(x,y) = cell_to_position(Cell::get_msx(), Cell::get_msy(), gid);
        geom_vec[2 * i]     = x;
        geom_vec[2 * i + 1] = y;
    }
}

int cpt_obj_size( void *data,
                  int num_gid_entries,
                  int num_lid_entries,
                  ZOLTAN_ID_PTR global_id,
                  ZOLTAN_ID_PTR local_id,
                  int *ierr) {
    ierr = ZOLTAN_OK;
    return sizeof(int) * 5 + sizeof(float);
}

void pack_particles(void *data,
                    int num_gid_entries,
                    int num_lid_entries,
                    ZOLTAN_ID_PTR global_id,
                    ZOLTAN_ID_PTR local_id,
                    int dest,
                    int size,
                    char *buf,
                    int *ierr) {
    auto all_mesh_data = (std::vector<Vehicle> *) data;
    memcpy(buf, &(all_mesh_data->operator[]((int)(*local_id))), sizeof(class std::vector<Vehicle>));
    all_mesh_data->operator[]((int)(*local_id)).gid = -1;
}

void unpack_particles ( void *data,
                        int num_gid_entries,
                        ZOLTAN_ID_PTR global_id,
                        int size,
                        char *buf,
                        int *ierr) {
    std::vector<Cell> *all_mesh_data = (std::vector<Vehicle> *) data;
    Vehicle v;
    memcpy(&v, buf, sizeof(int) * 5 + sizeof(float));
    all_mesh_data->push_back(v);
}

void post_migrate_particles (
        void *data,
        int num_gid_entries, int num_lid_entries, int num_import,
        ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
        int *import_procs, int num_export,
        ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids,
        int *export_procs, int *ierr) {
    auto all_mesh_data = (std::vector<Vehicle> *) data;
    size_t size = all_mesh_data->size();
    size_t i = 0;
    while(i < size) {
        if(all_mesh_data->operator[](i).gid == -1){
            std::iter_swap(all_mesh_data->begin() + i, all_mesh_data->end() - 1);
            all_mesh_data->pop_back();
            size--;
        } else {
            i++;
        }
    }
    size = all_mesh_data->size();
    for(size_t i = 0; i < size; i++){
        all_mesh_data->operator[](i).lid = i;
    }
}

Zoltan_Struct* zoltan_create_wrapper(bool automatic_migration, MPI_Comm comm, int num_global_part = -1, int part_on_me = -1) {
    std::string ngp = std::to_string(num_global_part);
    std::string pom = std::to_string(part_on_me);

    auto zz = Zoltan_Create(comm);

    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
    Zoltan_Set_Param(zz, "DETERMINISTIC", "1");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");

    if(num_global_part >= 1) Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", ngp.c_str());
    if(part_on_me >= 1)      Zoltan_Set_Param(zz, "NUM_LOCAL_PARTS",  pom.c_str());

    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");

    Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
    Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1");
    Zoltan_Set_Param(zz, "KEEP_CUTS", "1");

    if(automatic_migration)
        Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE");

    return zz;
}

Zoltan_Struct* zoltan_create_wrapper(bool automatic_migration = false) {
    return zoltan_create_wrapper(false, MPI_COMM_WORLD);
}

void zoltan_fn_init(Zoltan_Struct* zz, std::vector<Vehicle>* mesh_data, bool automatic_migration = false) {
    Zoltan_Set_Num_Obj_Fn(   zz, get_number_of_objects, mesh_data);
    Zoltan_Set_Obj_List_Fn(  zz, get_object_list,       mesh_data);
    Zoltan_Set_Num_Geom_Fn(  zz, get_num_geometry,      mesh_data);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list,     mesh_data);
    if(automatic_migration){
        Zoltan_Set_Obj_Size_Fn(zz, cpt_obj_size, mesh_data);
        Zoltan_Set_Pack_Obj_Fn(zz, pack_particles, mesh_data);
        Zoltan_Set_Unpack_Obj_Fn(zz, unpack_particles, mesh_data);
        Zoltan_Set_Post_Migrate_Fn(zz, post_migrate_particles, mesh_data);
    }
}

inline void zoltan_load_balance(std::vector<Vehicle>* mesh_data,
                                Zoltan_Struct* load_balancer,
                                bool automatic_migration = false,
                                bool do_migration = true) {

    /// ZOLTAN VARIABLES
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    /// END OF ZOLTAN VARIABLES

    automatic_migration = do_migration ? automatic_migration : false;

    zoltan_fn_init(load_balancer, mesh_data, automatic_migration);
    Zoltan_LB_Partition(load_balancer,      /* input (all remaining fields are output) */
                        &changes,           /* 1 if partitioning was changed, 0 otherwise */
                        &numGidEntries,     /* Number of integers used for a global ID */
                        &numLidEntries,     /* Number of integers used for a local ID */
                        &numImport,         /* Number of vertices to be sent to me */
                        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                        &importLocalGids,   /* Local IDs of vertices to be sent to me */
                        &importProcs,       /* Process rank for source of each incoming vertex */
                        &importToPart,      /* New partition for each incoming vertex */
                        &numExport,         /* Number of vertices I must send to other processes*/
                        &exportGlobalGids,  /* Global IDs of the vertices I must send */
                        &exportLocalGids,   /* Local IDs of the vertices I must send */
                        &exportProcs,       /* Process to which I send each of the vertices */
                        &exportToPart);     /* Partition to which each vertex will belong */

    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
}

#endif //SPEC_ZOLTAN_FN_HPP
