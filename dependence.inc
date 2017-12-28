driver_adjoint.o: constants.o
driver_adjoint.o: global_state.o
driver_adjoint.o: global.o
driver_adjoint.o: transit_to_solver.o
driver_adjoint.o: iteration_initialize.o
driver_adjoint.o: iteration_control.o
driver_adjoint.o: output_hdf5.o
driver_adjoint.o: output_visit.o

driver_adjoint_gpt.o: constants.o
driver_adjoint_gpt.o: global_state.o
driver_adjoint_gpt.o: global.o
driver_adjoint_gpt.o: transit_to_solver.o
driver_adjoint_gpt.o: iteration_initialize.o
driver_adjoint_gpt.o: iteration_control.o
driver_adjoint_gpt.o: reactivity.o
driver_adjoint_gpt.o: output_hdf5.o
driver_adjoint_gpt.o: output_visit.o

driver_perturb.o: constants.o
driver_perturb.o: global_state.o
driver_perturb.o: global.o
driver_perturb.o: driver_adjoint_gpt.o
driver_perturb.o: driver_steady.o
driver_perturb.o: transit_to_solver.o
driver_perturb.o: reactivity.o
driver_perturb.o: perturbation.o

driver_post_process.o: constants.o
driver_post_process.o: global_state.o
driver_post_process.o: global.o
driver_post_process.o: coefficient_iteration.o
driver_post_process.o: files_dirs.o
driver_post_process.o: output_main.o

driver_pre_process.o: constants.o
driver_pre_process.o: global_state.o
driver_pre_process.o: global.o
driver_pre_process.o: command_line.o
driver_pre_process.o: files_dirs.o
driver_pre_process.o: environment.o
driver_pre_process.o: input_driver.o

driver_steady.o: constants.o
driver_steady.o: global_state.o
driver_steady.o: global.o
driver_steady.o: stastics.o
driver_steady.o: transit_to_solver.o
driver_steady.o: iteration_initialize.o
driver_steady.o: iteration_control.o
driver_steady.o: feedback.o
driver_steady.o: time_advancing.o
driver_steady.o: output_timelist.o
driver_steady.o: output_hdf5.o
driver_steady.o: output_visit.o
driver_steady.o: test_interpolation.o

driver_transient_pcqs.o: constants.o
driver_transient_pcqs.o: global_state.o
driver_transient_pcqs.o: global.o
driver_transient_pcqs.o: driver_adjoint.o
driver_transient_pcqs.o: transit_to_solver.o
driver_transient_pcqs.o: transit_from_solver.o
driver_transient_pcqs.o: iteration_initialize.o
driver_transient_pcqs.o: iteration_control.o
driver_transient_pcqs.o: feedback.o
driver_transient_pcqs.o: time_advancing.o
driver_transient_pcqs.o: perturbation.o
driver_transient_pcqs.o: coefficient_pk_header.o
driver_transient_pcqs.o: process_pcqs.o
driver_transient_pcqs.o: process_pk.o
driver_transient_pcqs.o: output_timelist.o
driver_transient_pcqs.o: output_visit.o
driver_transient_pcqs.o: output_hdf5.o
driver_transient_pcqs.o: timestep_header.o

driver_transient_pk.o: constants.o
driver_transient_pk.o: global_state.o
driver_transient_pk.o: global.o
driver_transient_pk.o: transit_from_solver.o
driver_transient_pk.o: feedback.o
driver_transient_pk.o: time_advancing.o
driver_transient_pk.o: perturbation.o
driver_transient_pk.o: coefficient_pk_header.o
driver_transient_pk.o: process_pcqs.o
driver_transient_pk.o: process_pk.o
driver_transient_pk.o: output_timelist.o
driver_transient_pk.o: output_hdf5.o
driver_transient_pk.o: output_visit.o
driver_transient_pk.o: timestep_header.o

driver_transient_theta.o: constants.o
driver_transient_theta.o: global_state.o
driver_transient_theta.o: global.o
driver_transient_theta.o: transit_to_solver.o
driver_transient_theta.o: transit_from_solver.o
driver_transient_theta.o: iteration_initialize.o
driver_transient_theta.o: iteration_control.o
driver_transient_theta.o: feedback.o
driver_transient_theta.o: time_advancing.o
driver_transient_theta.o: perturbation.o
driver_transient_theta.o: process_pcqs.o
driver_transient_theta.o: output_timelist.o
driver_transient_theta.o: output_hdf5.o
driver_transient_theta.o: output_visit.o
driver_transient_theta.o: timestep_header.o

FRNK.o: global_state.o
FRNK.o: driver_pre_process.o
FRNK.o: driver_post_process.o
FRNK.o: driver_steady.o
FRNK.o: driver_adjoint.o
FRNK.o: driver_perturb.o
FRNK.o: driver_transient_theta.o
FRNK.o: driver_transient_pcqs.o
FRNK.o: driver_transient_pk.o
FRNK.o: th_pre_process.o
FRNK.o: th_post_process.o
FRNK.o: th_check_model.o

test_interpolation.o: constants.o
test_interpolation.o: global_state.o
test_interpolation.o: timer_header.o
test_interpolation.o: material_header.o
test_interpolation.o: global.o
test_interpolation.o: stastics.o
test_interpolation.o: input_xsec.o
test_interpolation.o: interpolation.o

command_line.o: constants.o

dict_header.o: constants.o

environment.o: constants.o

exception_header.o: constants.o

files_dirs.o: constants.o
files_dirs.o: string.o
files_dirs.o: exception_header.o

interpolation.o: constants.o
interpolation.o: exception_header.o

list_header.o: constants.o

search.o: constants.o
search.o: exception_header.o

set_header.o: constants.o
set_header.o: list_header.o

special_function.o: constants.o
special_function.o: exception_header.o

stastics.o: constants.o

string.o: constants.o

timer_header.o: constants.o

truncation.o: constants.o

vector_operation.o: constants.o
vector_operation.o: exception_header.o
vector_operation.o: set_header.o
vector_operation.o: string.o

adjoint_header.o: constants.o
adjoint_header.o: global_state.o
adjoint_header.o: geometry_header.o
adjoint_header.o: material_header.o

contain_header.o: constants.o
contain_header.o: global_state.o
contain_header.o: stastics.o
contain_header.o: geometry_header.o
contain_header.o: quadrature_header.o
contain_header.o: feedback_header.o

CRbank_header.o: constants.o
CRbank_header.o: exception_header.o
CRbank_header.o: state_header.o
CRbank_header.o: geometry_header.o
CRbank_header.o: material_header.o

decay_heat_header.o: constants.o
decay_heat_header.o: global_state.o
decay_heat_header.o: geometry_header.o
decay_heat_header.o: material_header.o

detector_header.o: constants.o
detector_header.o: timestep_header.o
detector_header.o: state_header.o
detector_header.o: geometry_header.o
detector_header.o: contain_header.o

feedback_header.o: constants.o
feedback_header.o: vector_operation.o
feedback_header.o: stastics.o
feedback_header.o: exception_header.o
feedback_header.o: state_header.o
feedback_header.o: geometry_header.o
feedback_header.o: link_header.o

geometry_header.o: constants.o
geometry_header.o: global_state.o

geomregular_header.o: constants.o
geomregular_header.o: string.o
geomregular_header.o: state_header.o
geomregular_header.o: geometry_header.o

global.o: constants.o
global.o: geometry_header.o
global.o: material_header.o
global.o: quadrature_header.o
global.o: sweeping_header.o
global.o: iteration_header.o
global.o: initvalue_header.o
global.o: adjoint_header.o
global.o: contain_header.o
global.o: CRbank_header.o
global.o: perturbation_header.o
global.o: worth_header.o
global.o: transient_header.o
global.o: link_header.o
global.o: feedback_header.o
global.o: LRAmodel_header.o
global.o: detector_header.o
global.o: precursor_solver_header.o
global.o: pksolver_header.o
global.o: pkreactivity_header.o

global_state.o: constants.o
global_state.o: state_header.o
global_state.o: timer_header.o
global_state.o: exception_header.o
global_state.o: timestep_header.o

initvalue_header.o: constants.o
initvalue_header.o: snap_interface.o
initvalue_header.o: iteration_header.o

iteration_header.o: constants.o
iteration_header.o: global_state.o
iteration_header.o: geometry_header.o
iteration_header.o: material_header.o
iteration_header.o: vector_operation.o

link_header.o: constants.o
link_header.o: Lilac_MODULE_Link.o

LRAmodel_header.o: constants.o
LRAmodel_header.o: exception_header.o
LRAmodel_header.o: state_header.o
LRAmodel_header.o: geometry_header.o
LRAmodel_header.o: material_header.o
LRAmodel_header.o: contain_header.o

material_header.o: constants.o
material_header.o: global_state.o
material_header.o: interpolation.o
material_header.o: vector_operation.o
material_header.o: exception_header.o
material_header.o: geometry_header.o

meshmap_header.o: constants.o

perturbation_header.o: constants.o
perturbation_header.o: global_state.o
perturbation_header.o: files_dirs.o
perturbation_header.o: material_header.o
perturbation_header.o: CRbank_header.o
perturbation_header.o: contain_header.o
perturbation_header.o: input_xsec.o

pkreactivity_header.o: constants.o

pksolver_header.o: constants.o
pksolver_header.o: lapack_interface.o
pksolver_header.o: state_header.o

precursor_solver_header.o: constants.o
precursor_solver_header.o: global_state.o
precursor_solver_header.o: timestep_header.o
precursor_solver_header.o: geometry_header.o
precursor_solver_header.o: material_header.o
precursor_solver_header.o: contain_header.o

quadrature_header.o: constants.o
quadrature_header.o: global_state.o
quadrature_header.o: state_header.o
quadrature_header.o: geometry_header.o

state_header.o: constants.o

sweeping_header.o: constants.o
sweeping_header.o: global_state.o
sweeping_header.o: geometry_header.o
sweeping_header.o: quadrature_header.o

timestep_header.o: constants.o
timestep_header.o: exception_header.o

transient_header.o: constants.o
transient_header.o: global_state.o

worth_header.o: constants.o

xesm_header.o: constants.o
xesm_header.o: global_state.o
xesm_header.o: material_header.o

hdf5_interface.o: constants.o

input_driver.o: constants.o
input_driver.o: global_state.o
input_driver.o: global.o
input_driver.o: th_global.o
input_driver.o: coefficient_iteration.o
input_driver.o: files_dirs.o
input_driver.o: input_keyword.o
input_driver.o: input_plain.o
input_driver.o: input_xml.o
input_driver.o: output_main.o
input_driver.o: output_runtime.o
input_driver.o: pkmodel_calculation.o

input_keyword.o: constants.o

input_plain.o: constants.o
input_plain.o: global_state.o
input_plain.o: geomregular_header.o
input_plain.o: global.o
input_plain.o: th_global.o
input_plain.o: string.o
input_plain.o: files_dirs.o
input_plain.o: input_keyword.o
input_plain.o: input_xsec.o
input_plain.o: output_main.o
input_plain.o: output_timelist.o
input_plain.o: output_visit.o

input_xml.o: constants.o
input_xml.o: global_state.o
input_xml.o: global.o
input_xml.o: th_global.o
input_xml.o: input_keyword.o
input_xml.o: xml_interface.o

input_xsec.o: constants.o
input_xsec.o: global_state.o
input_xsec.o: exception_header.o
input_xsec.o: interpolation.o
input_xsec.o: vector_operation.o
input_xsec.o: self_lagrange_interp.o
input_xsec.o: string.o
input_xsec.o: files_dirs.o
input_xsec.o: hdf5_interface.o
input_xsec.o: Lilac_MODULE_Link.o
input_xsec.o: material_header.o
input_xsec.o: CRbank_header.o
input_xsec.o: link_header.o
input_xsec.o: LRAmodel_header.o
input_xsec.o: feedback_header.o

lapack_interface.o: constants.o

Lilac_MODULE_Link.o: Lilac_MODULE_LinearSystem.o

output_hdf5.o: constants.o
output_hdf5.o: global_state.o
output_hdf5.o: global.o
output_hdf5.o: th_global.o
output_hdf5.o: hdf5_interface.o

output_main.o: constants.o
output_main.o: global_state.o
output_main.o: global.o

output_runtime.o: constants.o
output_runtime.o: global_state.o
output_runtime.o: global.o
output_runtime.o: coefficient_iteration.o

output_timelist.o: constants.o
output_timelist.o: global_state.o
output_timelist.o: contain_header.o

output_visit.o: constants.o
output_visit.o: global_state.o
output_visit.o: files_dirs.o
output_visit.o: global.o
output_visit.o: th_global.o

snap_interface.o: constants.o

xml_interface.o: constants.o

coefficient_pk_header.o: constants.o
coefficient_pk_header.o: global_state.o
coefficient_pk_header.o: global.o

feedback.o: constants.o
feedback.o: global_state.o
feedback.o: global.o
feedback.o: stastics.o
feedback.o: input_xsec.o
feedback.o: TH2NK_interface_self.o

perturbation.o: constants.o
perturbation.o: global_state.o
perturbation.o: global.o
perturbation.o: material_header.o

pkmodel_calculation.o: constants.o
pkmodel_calculation.o: global_state.o
pkmodel_calculation.o: global.o
pkmodel_calculation.o: th_global.o
pkmodel_calculation.o: stastics.o
pkmodel_calculation.o: timestep_header.o
pkmodel_calculation.o: output_timelist.o
pkmodel_calculation.o: driver_post_process.o
pkmodel_calculation.o: TH2NK_interface_self.o
pkmodel_calculation.o: th_pre_process.o
pkmodel_calculation.o: th_post_process.o

process_pcqs.o: constants.o
process_pcqs.o: global_state.o
process_pcqs.o: global.o
process_pcqs.o: timestep_header.o

process_pk.o: constants.o
process_pk.o: global_state.o
process_pk.o: lapack_interface.o
process_pk.o: global.o
process_pk.o: transient_header.o
process_pk.o: coefficient_pk_header.o

process_theta.o: constants.o
process_theta.o: global_state.o
process_theta.o: global.o

reactivity.o: constants.o
reactivity.o: global_state.o
reactivity.o: global.o
reactivity.o: contain_header.o

time_advancing.o: constants.o
time_advancing.o: global_state.o
time_advancing.o: global.o
time_advancing.o: transit_from_solver.o
time_advancing.o: timestep_header.o

transit_from_solver.o: constants.o
transit_from_solver.o: global_state.o
transit_from_solver.o: global.o

transit_to_solver.o: constants.o
transit_to_solver.o: global_state.o
transit_to_solver.o: global.o
transit_to_solver.o: process_theta.o

coefficient_header_accelarate.o: constants.o
coefficient_header_accelarate.o: global_state.o
coefficient_header_accelarate.o: iteration_header.o

coefficient_header_anisotropic.o: constants.o
coefficient_header_anisotropic.o: global_state.o
coefficient_header_anisotropic.o: quadrature_header.o

coefficient_header_transverse.o: constants.o
coefficient_header_transverse.o: global_state.o
coefficient_header_transverse.o: geometry_header.o
coefficient_header_transverse.o: material_header.o
coefficient_header_transverse.o: quadrature_header.o

coefficient_iteration.o: constants.o
coefficient_iteration.o: coefficient_header_transverse.o
coefficient_iteration.o: coefficient_header_anisotropic.o
coefficient_iteration.o: coefficient_header_accelarate.o

iteration_control.o: constants.o
iteration_control.o: global_state.o
iteration_control.o: iteration_header.o
iteration_control.o: global.o
iteration_control.o: coefficient_iteration.o
iteration_control.o: transit_from_solver.o
iteration_control.o: iteration_process.o
iteration_control.o: iteration_process_adjoint.o
iteration_control.o: output_main.o

iteration_initialize.o: constants.o
iteration_initialize.o: global_state.o
iteration_initialize.o: global.o

iteration_process.o: constants.o
iteration_process.o: global_state.o
iteration_process.o: truncation.o
iteration_process.o: global.o
iteration_process.o: coefficient_iteration.o
iteration_process.o: iteration_header.o
iteration_process.o: vector_operation.o

iteration_process_adjoint.o: constants.o
iteration_process_adjoint.o: global_state.o
iteration_process_adjoint.o: global.o
iteration_process_adjoint.o: coefficient_iteration.o
iteration_process_adjoint.o: iteration_header.o
iteration_process_adjoint.o: vector_operation.o

solver_1D.o: constants.o
solver_1D.o: global_state.o

abstract_property_header.o: constants.o

clad_steels_header.o: constants.o
clad_steels_header.o: exception_header.o
clad_steels_header.o: abstract_property_header.o

clad_Zr_header.o: constants.o
clad_Zr_header.o: exception_header.o
clad_Zr_header.o: abstract_property_header.o

coolant_HLM_header.o: constants.o
coolant_HLM_header.o: exception_header.o
coolant_HLM_header.o: abstract_property_header.o

coolant_water_parcs_header.o: constants.o
coolant_water_parcs_header.o: exception_header.o
coolant_water_parcs_header.o: abstract_property_header.o

coolant_water_refprop_header.o: constants.o
coolant_water_refprop_header.o: abstract_property_header.o

fuel_caramic_header.o: constants.o
fuel_caramic_header.o: exception_header.o
fuel_caramic_header.o: abstract_property_header.o

fuel_metallic_header.o: constants.o
fuel_metallic_header.o: exception_header.o
fuel_metallic_header.o: abstract_property_header.o

gap_gas_header.o: constants.o
gap_gas_header.o: exception_header.o
gap_gas_header.o: abstract_property_header.o

gth_geometry_header.o: constants.o

gth_parallel_kernel.o: constants.o
gth_parallel_kernel.o: lapack_interface.o
gth_parallel_kernel.o: stastics.o
gth_parallel_kernel.o: exception_header.o
gth_parallel_kernel.o: th_global.o
gth_parallel_kernel.o: abstract_property_header.o
gth_parallel_kernel.o: gth_geometry_header.o
gth_parallel_kernel.o: gth_thermal_header.o
gth_parallel_kernel.o: gth_power_header.o

gth_power_header.o: constants.o
gth_power_header.o: gth_geometry_header.o

gth_thermal_header.o: constants.o
gth_thermal_header.o: stastics.o
gth_thermal_header.o: abstract_property_header.o
gth_thermal_header.o: gth_geometry_header.o

TH2NK_interface_self.o: constants.o
TH2NK_interface_self.o: stastics.o
TH2NK_interface_self.o: th_global.o
TH2NK_interface_self.o: th_parallel_channel.o

th_check_model.o: constants.o
th_check_model.o: global.o
th_check_model.o: th_global.o
th_check_model.o: gth_parallel_kernel.o

th_check_property.o: constants.o
th_check_property.o: th_global.o
th_check_property.o: string.o
th_check_property.o: files_dirs.o

th_global.o: constants.o
th_global.o: coolant_water_parcs_header.o
th_global.o: coolant_water_refprop_header.o
th_global.o: coolant_HLM_header.o
th_global.o: clad_Zr_header.o
th_global.o: clad_steels_header.o
th_global.o: gap_gas_header.o
th_global.o: fuel_metallic_header.o
th_global.o: fuel_caramic_header.o
th_global.o: abstract_property_header.o
th_global.o: gth_geometry_header.o
th_global.o: gth_thermal_header.o
th_global.o: gth_power_header.o

th_output.o: constants.o
th_output.o: th_global.o

th_parallel_channel.o: constants.o
th_parallel_channel.o: global.o
th_parallel_channel.o: global_state.o
th_parallel_channel.o: th_global.o
th_parallel_channel.o: th_output.o
th_parallel_channel.o: gth_parallel_kernel.o

th_post_process.o: constants.o
th_post_process.o: global_state.o
th_post_process.o: files_dirs.o
th_post_process.o: th_global.o

th_pre_process.o: constants.o
th_pre_process.o: global.o
th_pre_process.o: th_global.o
th_pre_process.o: files_dirs.o

