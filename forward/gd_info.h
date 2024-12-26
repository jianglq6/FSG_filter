#ifndef GD_INFO_H
#define GD_INFO_H

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int ni;
  int nk;
  int nx;
  int nz;
  int ni1;
  int ni2;
  int nk1;
  int nk2;

  int lni;
  int lnk;
  int lnx;
  int lnz;
  int lni1;
  int lni2;
  int lnk1;
  int lnk2;

  int npoint_ghosts;
  int fdx_nghosts;
  int fdz_nghosts;

  // size of a single var
  //  the following two naming are same
  size_t siz_iz;
  size_t siz_icmp;

  // curvilinear coord name,
  char **index_name;
  
  //size_t siz_vars; // volume * num_of_vars, not easy for understand, may named with w3d and aux
} gdinfo_t;

/*************************************************
 * function prototype
 *************************************************/

int
gd_info_set_lebedev(
		gdinfo_t *const gdinfo,
        const int number_of_total_grid_points_x,
        const int number_of_total_grid_points_z,
              int abs_num_of_layers[][2],
        const int fdx_nghosts,
        const int fdz_nghosts,
        const int verbose);

int
gd_info_set(gdinfo_t *const gdinfo,
            const int number_of_total_grid_points_x,
            const int number_of_total_grid_points_z,
                  int abs_num_of_layers[][2],
            const int fdx_nghosts,
            const int fdz_nghosts,
            const int verbose);

int
gd_info_lindx_is_inner(int i, int k, gdinfo_t *gdinfo);

int
gd_info_pindx_is_inner(int i_phy, int k_phy, gdinfo_t *gdinfo);

int
gd_info_pindx_is_inner_i(int i_phy, gdinfo_t *gdinfo);

int
gd_info_pindx_is_inner_k(int k_phy, gdinfo_t *gdinfo);

int
gd_info_print(gdinfo_t *gdinfo);

#endif
