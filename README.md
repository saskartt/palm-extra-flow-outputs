# User code for additional diagnostic outputs for the PALM model system

This repo contains user code which implements various additional diagnostic outputs that can be useful for some studies. These are mainly helpful for computing both resolved and subgrid-scale Reynolds stress tensors, control volume studies and studies of surface drag force on buildings or other obstacles.

The implemented variables are:

* `pres_drag_norm_x*` - vertically integrated pressure drag force in x-direction computed from normalized perturbation pressure - kgm/s^2 - *only xy cross section*
* `pres_drag_norm_y*` - vertically integrated pressure drag force in y-direction computed from normalized perturbation pressure - kgm/s^2 - *only xy cross section*
* `wu_sgs`, `wv_sgs`, `uv_sgs`, `uw_sgs`, `vu_sgs`, `vw_sgs` - subgrid-scale components of the Reynolds stress tensor - m2/s2
* `utheta_product`, `vtheta_product` - products of the horizontal wind components and theta for total control volume flux - Km/s
* `uq_product`, `vq_product` - products of the horizontal wind components and q for total control volume flux - kgm/kgs
* `utheta_sqs`, `vtheta_sgs`, `wtheta_sgs` - subgrid-scale fluxes of sensible heat in three directions - W/m^2 (dynamic) or Km/s (kinematic)
* `uq_sqs`, `vq_sgs`, `wq_sgs` - subgrid-scale fluxes of moisture in three directions - W/m^2 (dynamic) or kgm/kgs (kinematic)

Both instantaneous and temporally averaged (`_av`) versions are implemented. Cross-sections `xy`, `xz` and `yz` are implemented as well.

Note: PALM implements both `pres_drag_x*` and `pres_drag_y*` natively. However, the results are incorrect with all-Neumann boundary conditions of perturbation pressure, as the reference is arbitrary. All-Neumann BCs for the perturbation pressure are used for nested domains. The outputs implemented here first normalize the perturbation pressure to zero mean before computing the pressure drag force.

## Usage
Drop `user_module.f90` into `USER_CODE` of your PALM job folder and add the desired inputs into `data_output_user` of the `&user_parameters` section of the `p3d` and `p3dr` namelists.

## Author
Sasu Karttunen \
<sasu.karttunen@helsinki.fi> \
University of Helsinki

Note that this is not a code that I would actively maintain and provide support for.
