from __future__ import print_function
from __future__ import absolute_import

import unittest

import numpy as np

from pyGDM2 import structures
from pyGDM2 import materials
from pyGDM2 import fields
from pyGDM2 import core
from pyGDM2 import linear
from pyGDM2 import tools




class TestClasses(unittest.TestCase):

    ## ---------- Setup test-simulations
    mesh1 = 'cube'
    mesh2 = 'hex'
    step = 30.0
    
    
    geometry1 = structures.sphere(step, R=5, mesh=mesh1)
    geometry1 = structures.center_struct(geometry1)
    geometry1.T[2] += 50
    n1, n2, n3 = 1.0, 1.0, 1.0
    spacing = 5000
    material1 = materials.silicon()
    struct1 = structures.struct(step, geometry1, material1, n1,n2, 
                    structures.get_normalization(mesh1), n3=n3, spacing=spacing)
    
    
    geometry2 = structures.sphere(step, R=5, mesh=mesh2)
    geometry2 = structures.center_struct(geometry2)
    geometry2.T[2] += 50
    n1, n2, n3 = 1.5, 1.2, 2.0
    spacing = 500
    material2 = materials.gold()
    struct2 = structures.struct(step, geometry2, material2, n1,n2, 
                    structures.get_normalization(mesh2), n3=n3, spacing=spacing)
    
    
    wavelengths = [700., 900.]
    field_generator = fields.planewave
    kwargs1 = dict(theta=[0.0, 90.0])
    efield1 = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs1)
    kwargs2 = dict(theta=np.linspace(0,90,20), kSign=[-1, 1])
    efield2 = fields.efield(field_generator, wavelengths=wavelengths, kwargs=kwargs2)
    
    
    ## ---------- simulation initialization
    simulations = []
    simulations.append(core.simulation(struct1, efield1))
    simulations.append(core.simulation(struct2, efield2))

        
# =============================================================================
#     expected results
# =============================================================================
    struct_norm = [[-0.00015514037795505151+4.821189765644286e-07j, -0.00015514037795505151+2.2684061585953231e-07j],
                   [-0.00015236224067758457+4.821189765644286e-07j, -0.00015236224067758457+2.2684061585953231e-07j]]
    struct_len_alpha = [[len(geometry1), len(geometry1)], [len(geometry2), len(geometry2)]]
    struct_alpha_first = [[29634.824 + 1j* 132.546, 26061.098 + 1j* 47.449795], [-27280.238 + 1j* 1617.3942, -51914.94 + 1j* 3036.9595]]
    
    efield_len_wavelengths = [2, 2]
    efield_len_kwargs = [1, 2]
    efield_len_kwargs_permutations = [2, 40]
# =============================================================================


    def test_0_simulation(self):
        print(" ----- Start core.simulation test -----")
        self.assertNotEqual(self.simulations[0], self.simulations[1])
        
        print("  sim.:", end='')
        for i_sim, sim in enumerate(self.simulations):
            print(" {},".format(i_sim), end='')
            self.assertTrue(sim.dtypef in [np.float32, np.float64])
            self.assertTrue(sim.dtypec in [np.complex64, np.complex128])
            
        print("\nFinished `simulation` test\n\n")


    def test_1_struct(self):
        print(" ----- Start structures.struct test -----")
        print("  sim.:", end='')
        for i_sim, sim in enumerate(self.simulations):
            print(" {},".format(i_sim), end='')
            for i_wl, wl in enumerate(sim.efield.wavelengths):
                norm = sim.struct.getNormalization(wl)
                self.assertAlmostEqual(norm, self.struct_norm[i_sim][i_wl], delta=1E-7)
                
                alpha = sim.struct.getPolarizability(wl)
                self.assertEqual(len(alpha), self.struct_len_alpha[i_sim][i_wl])
                self.assertAlmostEqual(alpha[0], self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                self.assertAlmostEqual(alpha.mean(), self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                
                alpha_tensor = sim.struct.getPolarizabilityTensor(wl)
                self.assertEqual(alpha_tensor[0].shape, (3,3))
                self.assertEqual(len(alpha_tensor), self.struct_len_alpha[i_sim][i_wl])
                self.assertAlmostEqual(alpha_tensor[0][0,0], self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                self.assertAlmostEqual(alpha_tensor[0][1,1], self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                self.assertAlmostEqual(alpha_tensor[0][2,2], self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                self.assertAlmostEqual(np.average(alpha_tensor, axis=0)[0,0], self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                self.assertAlmostEqual(np.average(alpha_tensor, axis=0)[1,1], self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                self.assertAlmostEqual(np.average(alpha_tensor, axis=0)[2,2], self.struct_alpha_first[i_sim][i_wl], delta=0.1)
                
#                print (norm, end='')
#                print (alpha[0].real, '+ 1j*', alpha[0].imag, end=', ')
#                print (np.array2string(alpha_tensor[0], separator=', '), end='\n')
                        
        print("\nFinished `struct` test\n\n")


    def test_2_efield(self):
        print(" ----- Start fields.efield test -----")
#        print("  sim.:", end='')

        for i_sim, sim in enumerate(self.simulations):
            print(" {},".format(i_sim), end='')
            self.assertEqual(len(sim.efield.wavelengths), self.efield_len_wavelengths[i_sim])
            self.assertEqual(len(sim.efield.kwargs), self.efield_len_kwargs[i_sim])
            self.assertEqual(len(sim.efield.kwargs_permutations), self.efield_len_kwargs_permutations[i_sim])

#            print(len(sim.efield.kwargs))
#            print(len(sim.efield.kwargs_permutations))
            
        print("\nFinished `efield` test\n\n")
        
        
    


class TestCore(unittest.TestCase):
    ## ---------- Setup test-simulations
    mesh = 'cube'
    step = 30.0
    
    geometry = structures.rect_wire(step, L=10, H=4, W=4, mesh=mesh)
    geometry = structures.center_struct(geometry)
    geometry.T[2] += 50
    
    
    
    n1, n2, n3 = 1.0, 1.0, 1.0
    spacing = 5000
    material1 = materials.dummy(2.5)
    struct1 = structures.struct(step, geometry, material1, n1,n2, 
                    structures.get_normalization(mesh), n3=n3, spacing=spacing)
    
    n1, n2, n3 = 1.5, 1.2, 2.0
    spacing = 500
    material2 = materials.silicon()
    struct2 = structures.struct(step, geometry, material2, n1,n2, 
                    structures.get_normalization(mesh), n3=n3, spacing=spacing)
    
    n1, n2, n3 = 2.5, 1.5, 1.0
    spacing = 400
    material3 = materials.gold()
    struct3 = structures.struct(step, geometry, material3, n1,n2, 
                    structures.get_normalization(mesh), n3=n3, spacing=spacing)
    
    n1, n2, n3 = 1.0, 1.5, 1.2
    spacing = 1500
    material4 = materials.silver()
    struct4 = structures.struct(step, geometry, material4, n1,n2, 
                    structures.get_normalization(mesh), n3=n3, spacing=spacing)
    
    
    wavelengths = [500.]                      # one single wavelength
    field_generator1 = fields.planewave
    kwargs1 = dict(theta=[0.0, 90.0])
    efield1 = fields.efield(field_generator1, wavelengths=wavelengths, kwargs=kwargs1)
    field_generator2 = fields.gaussian
    kwargs2 = dict(theta=[0.0, 90.0], spotsize=300)
    efield2 = fields.efield(field_generator2, wavelengths=wavelengths, kwargs=kwargs2)
    
    
    ## ---------- simulation initialization
    simulations = []
    for struct in [struct1, struct2, struct3, struct4]:
        simulations.append(core.simulation(struct, efield1))
    for struct in [struct1, struct2, struct3, struct4]:
        simulations.append(core.simulation(struct, efield2))
    for sim in simulations:
        sim.struct.getNormalization(sim.efield.wavelengths[0])
        sim.struct.getPolarizability(sim.efield.wavelengths[0])
    
    
# =============================================================================
#     expected results
# =============================================================================
    SBS_abs_sum = [6114.483 , 14756.186 , 3529.033 , 6873.047]
    K_abs_sum = [2069.2148 , 5252.8887 , 3431.2114 , 64333.88 , 2069.2148 , 5252.8887 , 3431.2114 , 64333.88]
    Escatter_abs_sum = [[ 201.93906 , 129.47525 ], [ 206.45422 , 105.12331 ], [ 109.123825 , 136.35785 ], [ 243.5806 , 323.61307 ], 
                        [ 225.66444 , 138.24681 ], [ 192.05823 , 102.62782 ], [ 92.782394 , 123.09842 ], [ 184.3924 , 305.2096 ]]
# =============================================================================


# =============================================================================
# CPU tests
# =============================================================================
    def test_0_get_sbs(self):
        print(" ----- Start core.get_side_by_side test -----")
        for i_sim, sim in enumerate(self.simulations[:4]):
            print(" {},".format(i_sim), end='')
            SBS = core.get_side_by_side(sim, sim.efield.wavelengths[0])
#            print(np.sum(np.abs(SBS.flatten())), ', ', end='')
            self.assertAlmostEqual(np.sum(np.abs(SBS.flatten())), self.SBS_abs_sum[i_sim], delta=.1)
            
            SBS = core.get_SBS_numba(sim, sim.efield.wavelengths[0])
            self.assertAlmostEqual(np.sum(np.abs(SBS.flatten())), self.SBS_abs_sum[i_sim], delta=.1)
            
        print("\nFinished get_sbs test\n\n")
        
    
    def test_1_get_K(self):
        print(" ----- Start core.get_general_propagator test -----", end='')
        methods = ["scipyinv", "numpyinv"]
        matrix_setup_methods = ["fortran", "numba"]
        for i_m, method in enumerate(methods):
            print("\n  testing '{}' inversion...".format(method), end='')
            for i_ms, matrix_setup in enumerate(matrix_setup_methods):
                print("\n    - with '{}' matrix setup method...".format(matrix_setup), end='')
                for i_sim, sim in enumerate(self.simulations):
                    print(" {},".format(i_sim), end='')
                    K = core.get_general_propagator(sim, wavelength=sim.efield.wavelengths[0],
                                                    method=method, verbose=0, 
                                                    matrix_setup=matrix_setup)
#                    print(np.sum(np.abs(K.flatten())), ', ', end='')
                    self.assertAlmostEqual(np.sum(np.abs(K.flatten())), self.K_abs_sum[i_sim], delta=.5)
        print("\nFinished get_general_propagator - CPU test\n\n")


    def test_2_scatter(self):
        print(" ----- Start core.scatter test -----", end='')
        methods = ["lu", "scipyinv"]
        matrix_setup_methods = ["fortran", "numba"]
        for i_m, method in enumerate(methods):
            print("\n  testing '{}' inversion...".format(method), end='')
            for i_ms, matrix_setup in enumerate(matrix_setup_methods):
                print("\n    - with '{}' matrix setup method...".format(matrix_setup))
                for i_sim, sim in enumerate(self.simulations):
                    if i_sim==0:
                        print("         testing plane wave illumination... ", end='')
                    if i_sim==4:
                        print("\n         testing gaussian illumination... ", end='')
                    print(" {},".format(i_sim), end='')
                    E = core.scatter(sim, method=method, matrix_setup=matrix_setup, verbose=0)
#                    print("[", np.sum(np.abs(E[0][1].flatten())), ', ', end='')
#                    print(np.sum(np.abs(E[1][1].flatten())), '], ', end='')
                    E_abs_sum = [np.sum(np.abs(E[0][1].flatten())), np.sum(np.abs(E[1][1].flatten()))]
                    self.assertAlmostEqual(E_abs_sum[0], self.Escatter_abs_sum[i_sim][0], delta=.002)
                    self.assertAlmostEqual(E_abs_sum[1], self.Escatter_abs_sum[i_sim][1], delta=.002)
        print("\nFinished scatter - CPU test\n\n")



# =============================================================================
# CUDA tests
# =============================================================================
    def test_GPU_0_get_sbs(self):
        print(" ----- Start core.get_side_by_side -  CUDA test -----")
        for i_sim, sim in enumerate(self.simulations[:4]):
            print(" {},".format(i_sim), end='')
            SBS = core.get_SBS_cuda(sim, sim.efield.wavelengths[0])
            self.assertAlmostEqual(np.sum(np.abs(SBS.flatten())), self.SBS_abs_sum[i_sim], delta=.1)
            
        print("\nFinished get_sbs CUDA test\n")


    def test_GPU_1_get_K(self):
        print(" ----- Start core.get_general_propagator -  CUDA test -----")
        for i_sim, sim in enumerate(self.simulations):
            print(" {},".format(i_sim), end='')
            K = core.get_general_propagator(sim, wavelength=sim.efield.wavelengths[0],
                                            method="cuda", 
                                            matrix_setup="fortran", verbose=0)
            K = K.get()
#            print(np.sum(np.abs(K.flatten())), ', ', end='')
            self.assertAlmostEqual(np.sum(np.abs(K.flatten())), self.K_abs_sum[i_sim], delta=.5)
                        
        print("\nFinished get_general_propagator - CUDA test\n\n")
    
    
    def test_GPU_2_scatter(self):
        print(" ----- Start core.scatter - CUDA test -----")
        for i_sim, sim in enumerate(self.simulations):
            if i_sim==0:
                print("         testing plane wave illumination... ", end='')
            if i_sim==4:
                print("\n         testing gaussian illumination... ", end='')
            print(" {},".format(i_sim), end='')
            E = core.scatter(sim, method="cuda", matrix_setup="cuda", verbose=0)
            E_abs_sum = [np.sum(np.abs(E[0][1].flatten())), np.sum(np.abs(E[1][1].flatten()))]
            self.assertAlmostEqual(E_abs_sum[0], self.Escatter_abs_sum[i_sim][0], delta=.002)
            self.assertAlmostEqual(E_abs_sum[1], self.Escatter_abs_sum[i_sim][1], delta=.002)
        print("\nFinished scatter - GPU test\n\n")
        




# =============================================================================
# run tests
# =============================================================================
if __name__ == '__main__':
    
#    ut_classes = TestClasses()
#    ut_classes.test_0_simulation()
#    ut_classes.test_1_struct()
#    ut_classes.test_2_efield()
    
#    ut_core = TestCore()
#    ut_core.test_0_get_sbs()
#    ut_core.test_1_get_K()
#    ut_core.test_2_scatter()
#    ut_core.test_GPU_0_get_sbs()
#    ut_core.test_GPU_1_get_K()
#    ut_core.test_GPU_2_scatter()
    
    
    ## all tests
    unittest.main()
    
    
    
    
    
    
    