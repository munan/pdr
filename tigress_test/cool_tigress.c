//#include "../copyright.h"
//TODO: uncomment the copyright include
/*==============================================================================
 * FILE: cool_tigress.c
 *
 * PURPOSE: heating and cooling function for TIGRESS simulations.
 *
 * REFERENCES: Gong, Ostriker and wolfire (2017)
 *
 * METHODS:
 * Estimate chemical abundances from gas and radiation parameters based on the
 * equilibrium chemistry in Gong, Ostriker and Wolfire (2017) network. Then
 * calculate:
 * Heating: cosmic-ray (CR) heating, photo-electric heating on dust (PE), UV
 * pumping of H2
 * Cooling: Ly-alpha, OI, C+, CI, CO, recombination of e- on PAHs
 * 
 * NOMENCLATURE:
 *
 *   Written by Munan Gong at 2018-10-15
 *
 * GLOBAL DEFAULT PARAMETERS
 * xCstd = 1.6e-4 (carbon abundance at solar metallicity)
 * xOstd = 3.2e-4 (oxygen abundance at solar metallicity)
 * xHetot = 0.1 (helium abundance)
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * get_abundances(nH, T, dvdr, Z, xi_CR, G_PE, G_CI, G_CO, G_H2,
 *                &x_e, &x_HI, &x_H2, &x_Cplus, &x_CI, &x_CO, &x_OI);
 * heating(x_e, x_HI, x_H2, nH, T, Z, xi_CR, G_PE, G_H2)
 * cooling(x_e, x_HI, x_H2, x_Cplus, x_CI, x_CO, x_OI, 
 *         nH, T, dvdr, Z, G_PE)
 *
 * Notes:
 * - Input parameters nH, T, dvdr are in CGS units, radiation field G_i in
 * units of Draine 1978 radiation field, xi_CR in [s^{-1}H^{-1}] (per second per
 * H atom).
 * - Output heating and cooling rates in [ergs s^{-1} cm^3] units (note that
 * the private heating and cooling functions such as heatingPE, coolingLya,
 * etc. are in units of [ergs s^{-1} H^{-1}], differs from the public heating()
 * and cooling() functions by nH).
 * - Z = Z_d = Z_g is the same for gas and dust metallicity
 *
 * CONTAINS PRIVATE FUNCTIONS:
 * chemical abundances:
 * fH2(nH, T, Z_d, xi_CR, G_H2)
 * fCplus(x_e, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)
 * fHplus(x_e, x_H2, x_Cplus, nH, T, Z_d, xi_CR, G_PE)
 * fe_e(x_e, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)
 * fe(x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)
 * fCO(x_H2, x_Cplus, nH, Z_d, Z_g, xi_CR, G_CO)
 *
 * Heating:
 * heatingPE(x_e, nH, T, Z_d, G_PE)
 * heatingCR(x_e, x_HI, x_H2, nH, xi_CR)
 * heatingH2pump(x_HI, x_H2, nH, T, G_H2)
 * 
 * Cooling:
 * coolingLya(x_e, x_HI, nH, T)
 * coolingOI(x_e, x_OI, x_HI, x_H2, nH, T)
 * coolingCII(x_e, x_Cplus, x_HI, x_H2, nH, T)
 * coolingCI(x_e, x_CI, x_HI, x_H2, nH, T)
 * coolingCO(x_e, x_CO, x_HI, x_H2, nH, T, dvdr)
 * coolingRec(x_e, nH, T, Z_d, G_PE)
 * coolingHot(T, Z_g)
 *
 * Notes:
 * The heating and cooling rates above are in [ergs s^{-1} H^{-1}] units. 
 * In order to get the heating and cooling rates for TIGRESS in
 * [ergs s^{-1} cm^{-3}] units, they are multiplied by density nH in functions
 * heating() and cooling().
 *
 * helper functions for cooling:
 * CII_rec_rate_()
 * q10CII_()
 * cooling2Level_()
 * cooling3Level_()
 * linearInterpIndex_()
 * linearInterp_()
 * LP1Di_()
 * LP2Di_()
 *
 *============================================================================*/

#include "cool_tigress.h" //TODO: merge header into prototypes.h
//TODO: include ../defs.h which has definition of MIN and MAX, then delete the
//two definitions below
#define MAX(x, y) (((x) > (y)) ? (x) : (y)) 
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/*elemental abundances*/
static const Real xCstd=1.6e-4, xOstd=3.2e-4, xHetot=0.1;
/*physical constants*/
static const Real eV_ = 1.602e-12;
static const Real kb_ = 1.381e-16;
/*ortho to para ratio of H2*/
static const Real fo_ = 0.75;
static const Real fp_ = 0.25;
/*photo-electric heating*/
static const Real CPE_[7] = {5.22, 2.25, 0.04996, 0.00430, 0.147, 0.431,0.692};
static const Real DPE_[5] = {0.4535, 2.234, -6.266, 1.442, 0.05089};
/*----C+, 2 level system---*/
static const Real A10CII_ = 2.3e-6;
static const Real E10CII_ = 1.26e-14;
static const Real g0CII_ = 2.;
static const Real g1CII_ = 4.;
/*----HI, 2 level system---*/
static const Real A10HI_ = 6.265e8;
static const Real E10HI_ = 1.634e-11;
static const Real g0HI_ = 1.;
static const Real g1HI_ = 3.;
/*----CI, 3 level system---*/
static const Real g0CI_ = 1;
static const Real g1CI_ = 3;
static const Real g2CI_ = 5;
static const Real A10CI_ = 7.880e-08;
static const Real A20CI_ = 1.810e-14;
static const Real A21CI_ = 2.650e-07;
static const Real E10CI_ = 3.261e-15;
static const Real E20CI_ = 8.624e-15;
static const Real E21CI_ = 5.363e-15;
/*----OI, 3 level system---*/
static const Real g0OI_ = 5;
static const Real g1OI_ = 3;
static const Real g2OI_ = 1;
static const Real A10OI_ = 8.910e-05;
static const Real A20OI_ = 1.340e-10;
static const Real A21OI_ = 1.750e-05;
static const Real E10OI_ = 3.144e-14;
static const Real E20OI_ = 4.509e-14;
static const Real E21OI_ = 1.365e-14;

/*-----CO cooling table data, from Omukai+2010-----*/
static const int lenTCO_ = 11;
static const int lenNeffCO_ = 11;
static const Real TCO_[lenTCO_] = {10,	20,	30,	50,	80,	100,
                                   300,	600,	1000,	1500,	2000};
static const Real NeffCO_[lenNeffCO_] = {14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
                                         17.0, 17.5, 18.0, 18.5, 19.0};
static const Real L0CO_[lenTCO_] = {24.77,	24.38, 24.21,	24.03, 23.89, 23.82,
  23.34238089,  22.99832519,  22.75384686,  22.56640625, 22.43740866};
static const Real LLTECO_[lenNeffCO_*lenTCO_] = {
21.08, 20.35, 19.94, 19.45, 19.01, 18.80, 17.81, 17.23, 16.86, 16.66, 16.55,
21.09, 20.35, 19.95, 19.45, 19.01, 18.80, 17.81, 17.23, 16.86, 16.66, 16.55,
21.11, 20.37, 19.96, 19.46, 19.01, 18.80, 17.81, 17.23, 16.86, 16.66, 16.55,
21.18, 20.40, 19.98, 19.47, 19.02, 18.81, 17.82, 17.23, 16.87, 16.66, 16.55,
21.37, 20.51, 20.05, 19.52, 19.05, 18.83, 17.82, 17.23, 16.87, 16.66, 16.55,
21.67, 20.73, 20.23, 19.64, 19.13, 18.90, 17.85, 17.25, 16.88, 16.67, 16.56,
22.04, 21.05, 20.52, 19.87, 19.32, 19.06, 17.92, 17.28, 16.90, 16.69, 16.58,
22.44, 21.42, 20.86, 20.19, 19.60, 19.33, 18.08, 17.38, 16.97, 16.75, 16.63,
22.87, 21.82, 21.24, 20.55, 19.95, 19.66, 18.34, 17.59, 17.15, 16.91, 16.78,
23.30, 22.23, 21.65, 20.94, 20.32, 20.03, 18.67, 17.89, 17.48, 17.26, 17.12,
23.76, 22.66, 22.06, 21.35, 20.71, 20.42, 19.03, 18.26, 17.93, 17.74, 17.61
};
static const Real nhalfCO_[lenNeffCO_*lenTCO_] = {
  3.29, 3.49 ,3.67  ,3.97,  4.30, 4.46, 5.17, 5.47, 5.53, 5.30, 4.70, 
  3.27, 3.48 ,3.66  ,3.96,  4.30, 4.45, 5.16, 5.47, 5.53, 5.30, 4.70, 
  3.22, 3.45 ,3.64  ,3.94,  4.29, 4.45, 5.16, 5.47, 5.53, 5.30, 4.70, 
  3.07, 3.34 ,3.56  ,3.89,  4.26, 4.42, 5.15, 5.46, 5.52, 5.30, 4.70, 
  2.72, 3.09 ,3.35  ,3.74,  4.16, 4.34, 5.13, 5.45, 5.51, 5.29, 4.68, 
  2.24, 2.65 ,2.95  ,3.42,  3.92, 4.14, 5.06, 5.41, 5.48, 5.26, 4.64, 
  1.74, 2.15 ,2.47  ,2.95,  3.49, 3.74, 4.86, 5.30, 5.39, 5.17, 4.53, 
  1.24, 1.65 ,1.97  ,2.45,  3.00, 3.25, 4.47, 5.02, 5.16, 4.94, 4.27, 
 0.742, 1.15 ,1.47  ,1.95,  2.50, 2.75, 3.98, 4.57, 4.73, 4.52, 3.84, 
 0.242, 0.652,0.966 ,1.45,  2.00, 2.25, 3.48, 4.07, 4.24, 4.03, 3.35, 
-0.258, 0.152,0.466 ,0.95,	 1.50, 1.75, 2.98, 3.57, 3.74, 3.53, 2.85
};
static const Real alphaCO_[lenNeffCO_*lenTCO_] = {
0.439, 0.409, 0.392, 0.370, 0.361, 0.357, 0.385, 0.437, 0.428, 0.354, 0.322,	
0.436, 0.407, 0.391, 0.368, 0.359, 0.356, 0.385, 0.437, 0.427, 0.354, 0.322,	
0.428, 0.401, 0.385, 0.364, 0.356, 0.352, 0.383, 0.436, 0.427, 0.352, 0.320,	
0.416, 0.388, 0.373, 0.353, 0.347, 0.345, 0.380, 0.434, 0.425, 0.349, 0.316,	
0.416, 0.378, 0.360, 0.338, 0.332, 0.330, 0.371, 0.429, 0.421, 0.341, 0.307,	
0.450, 0.396, 0.367, 0.334, 0.322, 0.317, 0.355, 0.419, 0.414, 0.329, 0.292,	
0.492, 0.435, 0.403, 0.362, 0.339, 0.329, 0.343, 0.406, 0.401, 0.317, 0.276,	
0.529, 0.473, 0.441, 0.404, 0.381, 0.370, 0.362, 0.410, 0.392, 0.316, 0.272,	
0.555, 0.503, 0.473, 0.440, 0.423, 0.414, 0.418, 0.446, 0.404, 0.335, 0.289,	
0.582, 0.528, 0.499, 0.469, 0.457, 0.451, 0.470, 0.487, 0.432, 0.364, 0.310,	
0.596, 0.546, 0.519, 0.492, 0.483, 0.479, 0.510, 0.516, 0.448, 0.372, 0.313
};

/*-----Cooling for hot gas, Wiersma, Schaye, and Smith (2008)-----*/
static const int len_hot_ = 262;
static const Real log10T_hot_[len_hot_] = {
 3.79022181, 3.81011149, 3.82999809, 3.84989215, 3.86978304, 3.88967699,
 3.90956673, 3.9294598 , 3.94935097, 3.9692388 , 3.98913377, 4.00902574,
 4.02889645, 4.04879127, 4.06870522, 4.08859673, 4.10846354, 4.12836695,
 4.14826323, 4.16814378, 4.18805621, 4.20793044, 4.22783531, 4.24772783,
 4.26761753, 4.28751061, 4.30738906, 4.32727718, 4.34717384, 4.36705759,
 4.38696244, 4.40684663, 4.42673892, 4.44663035, 4.46652657, 4.48641631,
 4.50630204, 4.52619707, 4.54608592, 4.56597762, 4.58586663, 4.60575743,
 4.62565193, 4.64554998, 4.66544027, 4.68532963, 4.70522206, 4.72511088,
 4.74500449, 4.76489314, 4.78478128, 4.80467771, 4.82456837, 4.84445854,
 4.86435086, 4.88424011, 4.90413104, 4.92402587, 4.94391476, 4.96380671,
 4.98369844, 5.00358977, 5.02349938, 5.04336228, 5.06325828, 5.08314414,
 5.10305075, 5.12293637, 5.14282729, 5.16271373, 5.18261435, 5.20251556,
 5.22240429, 5.2422929 , 5.26216621, 5.28207805, 5.30196273, 5.32184688,
 5.34175098, 5.36163341, 5.38153022, 5.40141777, 5.42130772, 5.44119223,
 5.46109316, 5.48098373, 5.50086742, 5.52075856, 5.54065476, 5.56054026,
 5.58043455, 5.60033023, 5.62021936, 5.64011357, 5.66000172, 5.67989102,
 5.69978571, 5.71967925, 5.73956443, 5.75945631, 5.77935114, 5.79924403,
 5.81913577, 5.83902529, 5.85891599, 5.87880893, 5.89869776, 5.9185912 ,
 5.9384797 , 5.95837273, 5.97826271, 5.99815468, 6.0180344 , 6.03794416,
 6.05781819, 6.07773118, 6.09760433, 6.1175033 , 6.13738576, 6.1572754 ,
 6.17719008, 6.19705991, 6.21695721, 6.23683945, 6.25674179, 6.27662262,
 6.29653357, 6.31641071, 6.3362996 , 6.35619801, 6.37608399, 6.39597255,
 6.41587441, 6.43576476, 6.45565175, 6.47554044, 6.49543332, 6.51533064,
 6.53521814, 6.55510655, 6.575003  , 6.59488955, 6.61478125, 6.63467875,
 6.65457118, 6.67445696, 6.69435069, 6.71424591, 6.73413548, 6.75402708,
 6.77391803, 6.79380435, 6.81370105, 6.83359329, 6.85347918, 6.87337292,
 6.89326229, 6.91315662, 6.9330467 , 6.95293767, 6.97282744, 6.99272137,
 7.01262635, 7.03249788, 7.0523861 , 7.07228667, 7.09219412, 7.1120685 ,
 7.13197135, 7.15185996, 7.1717557 , 7.1916466 , 7.211521  , 7.23141861,
 7.25129746, 7.27119084, 7.29108011, 7.31099053, 7.33088029, 7.35077118,
 7.3706611 , 7.39054654, 7.41043979, 7.4303331 , 7.45021831, 7.47011635,
 7.49000064, 7.50990113, 7.52978955, 7.54967749, 7.56957287, 7.58945809,
 7.60934892, 7.62924645, 7.64913032, 7.6690285 , 7.68891791, 7.70881167,
 7.7287027 , 7.74859111, 7.7684827 , 7.78837345, 7.808265  , 7.82815701,
 7.84804741, 7.86793865, 7.88783102, 7.90772299, 7.92761157, 7.94750234,
 7.96739351, 7.98728631, 8.00719282, 8.02706406, 8.04696315, 8.06684751,
 8.08675123, 8.10663279, 8.12652103, 8.14640714, 8.16631168, 8.18619325,
 8.20609694, 8.22598088, 8.24588265, 8.26576092, 8.28564731, 8.30554482,
 8.32543356, 8.34533451, 8.36522583, 8.38510556, 8.40500467, 8.42489796,
 8.44479401, 8.46468325, 8.4845703 , 8.50445727, 8.52435717, 8.54424173,
 8.56413322, 8.58402575, 8.60391264, 8.62380731, 8.64369936, 8.66358786,
 8.68347932, 8.70337737, 8.72326681, 8.74315686, 8.76304595, 8.7829382 ,
 8.80282844, 8.82272367, 8.84261548, 8.86250705, 8.8823936 , 8.9022858 ,
 8.9221803 , 8.94206761, 8.96196185, 8.98185031
};

static const Real log10L_hot_HHe_[len_hot_] = {
 -28.22599427, -28.17063189, -28.11013268, -28.03914409, -27.9486545 , 
 -27.82684288, -27.65908021, -27.43191567, -27.1340969 , -26.76376514,
 -26.3302091 , -25.8562048 , -25.36138097, -24.85942911, -24.35726846,
 -23.87127772, -23.42018212, -23.01582386, -22.66669382, -22.38175312,
 -22.16586553, -22.01895648, -21.93081307, -21.8907928 , -21.88352509,
 -21.89973287, -21.92903908, -21.96734095, -22.00878279, -22.05216942,
 -22.09347102, -22.13220264, -22.16769551, -22.20042513, -22.23081475,
 -22.26233037, -22.29519767, -22.32956859, -22.36402393, -22.39970247,
 -22.43415218, -22.46636773, -22.49284604, -22.51575712, -22.53035468,
 -22.53345859, -22.52038151, -22.49180707, -22.4454745 , -22.38329348,
 -22.30883891, -22.22995107, -22.15230398, -22.08098993, -22.02283906,
 -21.98163242, -21.95809436, -21.95140284, -21.95916041, -21.97959728,
 -22.00794159, -22.04190459, -22.07952187, -22.1194951 , -22.16067282,
 -22.20229631, -22.2437748 , -22.28484938, -22.32521232, -22.3649085 ,
 -22.40329297, -22.44087173, -22.47732108, -22.51243744, -22.54660521,
 -22.57968112, -22.61163233, -22.64250348, -22.67227335, -22.701016  ,
 -22.72871613, -22.7554013 , -22.78107081, -22.80579102, -22.82953292,
 -22.87998575, -22.90305524, -22.92492774, -22.94577009, -22.96569247,
 -22.98464024, -23.00269227, -23.01981931, -23.03604698, -23.05140612,
 -23.06591079, -23.07966271, -23.09449555, -23.10836778, -23.12131164,
 -23.13334057, -23.14449902, -23.15480891, -23.16430309, -23.17296045,
 -23.18087082, -23.18805654, -23.19454672, -23.2003838 , -23.20559561,
 -23.21022193, -23.21427876, -23.21686855, -23.2188049 , -23.22020877,
 -23.2211183 , -23.22155932, -23.22154485, -23.22112553, -23.22030973,
 -23.21913566, -23.21761338, -23.2157751 , -23.21363905, -23.21123082,
 -23.208562  , -23.2056514 , -23.20252471, -23.19918673, -23.19565608,
 -23.19195791, -23.18809002, -23.18407708, -23.17992339, -23.17564634,
 -23.17125635, -23.16675724, -23.16089922, -23.15454846, -23.14809375,
 -23.14154101, -23.13489603, -23.1281703 , -23.12136907, -23.11450308,
 -23.10757702, -23.10060087, -23.09357875, -23.08651464, -23.07942278,
 -23.07230122, -23.06515324, -23.05799195, -23.05081987, -23.04363464,
 -23.03644827, -23.02926271, -23.02207529, -23.0149015 , -23.00773375,
 -23.000578  , -22.99344839, -22.98518555, -22.97629496, -22.96738124,
 -22.95848889, -22.94958157, -22.9407396 , -22.93188838, -22.92303205,
 -22.91421011, -22.90542407, -22.89664106, -22.88789795, -22.87919582,
 -22.87053569, -22.86188685, -22.85325199, -22.84469433, -22.83615312,
 -22.82765983, -22.81921487, -22.8108186 , -22.80247134, -22.7941463 ,
 -22.78589826, -22.77767379, -22.76865786, -22.75877737, -22.74897046,
 -22.73916569, -22.72941432, -22.71971763, -22.71003221, -22.70040555,
 -22.69081713, -22.68126903, -22.67176315, -22.66232115, -22.65290427,
 -22.64355344, -22.63425002, -22.62499518, -22.61579   , -22.60663541,
 -22.59753227, -22.58849814, -22.57951615, -22.57058689, -22.56171087,
 -22.55290402, -22.54415057, -22.53485518, -22.52457581, -22.51432115,
 -22.50410892, -22.49394167, -22.48382173, -22.47373828, -22.46369413,
 -22.45370416, -22.44374576, -22.43384545, -22.42398132, -22.41417845,
 -22.40440459, -22.39469495, -22.385029  , -22.37540854, -22.3658453 ,
 -22.35632037, -22.34685505, -22.33743103, -22.32805926, -22.31874067,
 -22.30946725, -22.3002403 , -22.29071307, -22.28000617, -22.26932268,
 -22.25865928, -22.24802074, -22.23740395, -22.2268135 , -22.2162393 ,
 -22.20568626, -22.19515209, -22.1846348 , -22.17413267, -22.16364423,
 -22.15316825, -22.1427037 , -22.13224975, -22.12179427, -22.11134321,
 -22.10089114, -22.09043862, -22.07997064, -22.06949414, -22.05900068,
 -22.04847792, -22.03792914
};

static const Real log10L_hot_metal_[len_hot_] = {
 -25.575347  , -25.56504806, -25.55453591, -25.54322416, -25.53039887,
 -25.51467615, -25.4953117 , -25.47114755, -25.44456696, -25.41715289,
 -25.3837451 , -25.32151804, -25.19078382, -24.97444838, -24.68552182,
 -24.36203021, -24.03941032, -23.74550026, -23.47395811, -23.23994713,
 -23.04767745, -22.90097578, -22.79490395, -22.72402101, -22.67458697,
 -22.63815888, -22.60488338, -22.57265998, -22.53809648, -22.50258654,
 -22.46442741, -22.42433885, -22.38251738, -22.3392009 , -22.29486357,
 -22.25010143, -22.20538649, -22.160981  , -22.11801083, -22.07482203,
 -22.03204974, -21.99003947, -21.94861592, -21.90805181, -21.86828512,
 -21.82932766, -21.79115571, -21.75384623, -21.71735523, -21.68166532,
 -21.64677642, -21.61280601, -21.57973062, -21.54643706, -21.51524475,
 -21.48610351, -21.45918269, -21.43437647, -21.41140339, -21.39126008,
 -21.37300304, -21.35820871, -21.34827629, -21.34382511, -21.34405522,
 -21.34676818, -21.34921667, -21.3492749 , -21.34603536, -21.33974182,
 -21.33166068, -21.32263033, -21.31447322, -21.30783192, -21.30256141,
 -21.29853564, -21.29463242, -21.29045021, -21.2862911 , -21.28418879,
 -21.28644226, -21.29633898, -21.31731352, -21.35150338, -21.39905976,
 -21.4558327 , -21.5153301 , -21.57071617, -21.61769472, -21.65415596,
 -21.68041513, -21.69808061, -21.70925328, -21.71615624, -21.72033306,
 -21.72367609, -21.72764776, -21.73367421, -21.74241741, -21.75483293,
 -21.76942336, -21.78457395, -21.7977297 , -21.80723746, -21.81284605,
 -21.81423127, -21.81253568, -21.80924822, -21.80559667, -21.80274702,
 -21.80128794, -21.80222339, -21.80556892, -21.81107152, -21.81841364,
 -21.82728034, -21.83752521, -21.84921287, -21.86264589, -21.87804642,
 -21.89650416, -21.91850867, -21.94469814, -21.9753912 , -22.01048789,
 -22.04922636, -22.09036375, -22.13220264, -22.1738215 , -22.21401566,
 -22.25232412, -22.28885948, -22.3240934 , -22.3587044 , -22.3931459 ,
 -22.42726777, -22.46120013, -22.49403648, -22.52511387, -22.55385132,
 -22.57986265, -22.60292945, -22.62296709, -22.639937  , -22.6539993 ,
 -22.66526485, -22.67399203, -22.68051917, -22.68510098, -22.68790731,
 -22.68915805, -22.68896703, -22.68727272, -22.68486968, -22.68191579,
 -22.67875379, -22.67592342, -22.6737256 , -22.6724776 , -22.67251846,
 -22.67415607, -22.67769799, -22.6833588 , -22.69143559, -22.70215198,
 -22.71577236, -22.73254666, -22.75249318, -22.77564825, -22.801893  ,
 -22.83088422, -22.86217134, -22.89513748, -22.92903908, -22.96329128,
 -22.99718622, -23.02949506, -23.06088528, -23.09062586, -23.11856398,
 -23.14465654, -23.16901378, -23.19152569, -23.21227039, -23.23132491,
 -23.2487825 , -23.26479246, -23.2792455 , -23.29210635, -23.30366983,
 -23.31401656, -23.32323325, -23.33139978, -23.33856595, -23.34479431,
 -23.35013963, -23.354676  , -23.35845649, -23.36152073, -23.36390345,
 -23.3656739 , -23.36687508, -23.36756296, -23.36778574, -23.3675832 ,
 -23.36701661, -23.36610756, -23.36497894, -23.36377296, -23.36241021,
 -23.3609221 , -23.35938939, -23.35744645, -23.35481361, -23.35211846,
 -23.34946907, -23.34676818, -23.34406481, -23.34134956, -23.33865119,
 -23.33593182, -23.33319196, -23.33042281, -23.32761579, -23.32476257,
 -23.32185507, -23.31886735, -23.31581016, -23.31266727, -23.30943183,
 -23.30608847, -23.3026573 , -23.29911464, -23.29546342, -23.2917066 ,
 -23.28784716, -23.28387975, -23.27947687, -23.27379458, -23.2679844 ,
 -23.26204443, -23.25598098, -23.24980799, -23.24351604, -23.23712646,
 -23.23063006, -23.22408387, -23.21741989, -23.21068728, -23.20387685,
 -23.19700064, -23.19005666, -23.18306306, -23.17601744, -23.16892407,
 -23.16179969, -23.15462906, -23.14743462, -23.14020745, -23.1329511 ,
 -23.12567477, -23.11837572
};


/*----------------------------------------------------------------------------*/
/* PRIVATE FUCNTIONS                                                          */
/*----------------------------------------------------------------------------*/
//Chemistry-------------------------------------------------------------------
//H2 abundance consider CR destruction and FUV
static Real fH2(const Real nH, const Real T, const Real Z_d,
                const Real xi_CR, const Real G_H2);
//C+ abundance, depending on abundances of e- and H2
static Real fCplus(const Real x_e, const Real x_H2, const Real nH, const Real T,
                   const Real Z_d, const Real Z_g, const Real xi_CR, 
                   const Real G_PE, const Real G_CI);
//H+ abundance, depending on abundances of C+, e- and H2
static Real fHplus(const Real x_e, const Real x_Cplus, const Real x_H2,
                   const Real nH, const Real T, const Real Z_d, const Real xi_CR,
                   const Real G_PE);
//ion abundances: ion = C+ + H+, depending on e- and H2 abundances.
static Real fions(const Real x_e, const Real x_H2, const Real nH, const Real T,
                  const Real Z_d, const Real Z_g, const Real xi_CR, 
                  const Real G_PE, const Real G_CI);
//e- abundances, solved from iteration, assume e- balances from C+ and H+
static Real fe(const Real x_H2, const Real nH, const Real T,
               const Real Z_d, const Real Z_g, const Real xi_CR, 
               const Real G_PE, const Real G_CI);
//CO abundance. Fit from Gong, Ostriker and Wolfire 2017
//Note: This is only tested in the case Z_d = Z_g
static Real fCO(const Real x_H2, const Real x_Cplus, const Real nH, 
                const Real Z_d, const Real Z_g, const Real xi_CR, const Real G_CO);

/* TODO: move this back to private functions
//heating---------------------------------------------------------------------
//cosmic ray heating
static Real heaingCR(const Real x_e, const Real x_HI, const Real x_H2,  
                     const Real nH, const Real xi_CR);
//photo electric heating on dust
static Real heatingPE(const Real x_e, const Real nH, const Real T, 
                      const Real Z_d, const Real G_PE);
//UV-pumping of H2
static Real heatingH2pump(const Real x_HI, const Real x_H2, const Real nH,
                          const Real T, const Real G_H2);

//cooling---------------------------------------------------------------------
//HI Lyman alpha cooling
static Real coolingLya(const Real x_e, const Real x_HI, const Real nH,
                       const Real T);
//OI cooling
static Real coolingOI(const Real x_e, const Real x_OI, const Real x_HI,
                      const Real x_H2, const Real nH, const Real T);
//C+ cooling
static Real coolingCII(const Real x_e, const Real x_Cplus, const Real x_HI,
                       const Real x_H2, const Real nH, const Real T);
//CI cooling
static Real coolingCI(const Real x_e, const Real x_CI, const Real x_HI,
                      const Real x_H2, const Real nH, const Real T);
//CO rotational line cooling, dvdr in cgs units
static Real coolingCO(const Real x_e, const Real x_CO, const Real x_HI,
                      const Real x_H2, const Real nH, const Real T, 
                      const Real dvdr);
//cooling by recombination of e on PAHs
static Real coolingRec(const Real x_e, const Real nH, const Real T, 
                       const Real Z_d, const Real G_PE);
*/

//cooling of hot gas, Wiersma, Schaye, and Smith (2008)
static Real coolingHot(const Real T, const Real Z_g);

//helper functions for cooling-------------------------------------------------
static Real CII_rec_rate_(const Real temp);
static Real q10CII_(const Real nHI, const Real nH2, const Real ne, const Real T);
static Real cooling2Level_(const Real q01, const Real q10, const Real A10,
                           const Real E10, const Real xs);

static Real cooling3Level_(const Real q01, const Real q10, const Real q02,
                           const Real q20, const Real q12, const Real q21,
													 const Real A10, const Real A20, const Real A21,
                           const Real E10, const Real E20, const Real E21,
													 const Real xs);
static int linearInterpIndex_(const int len, const Real xarr[], const Real x);
static Real linearInterp_(const Real x0, const Real x1, const Real y0,
                          const Real y1, const Real x);
static Real LP1Di_(const Real *xarr, const Real *data, const int ix,
                   const Real x);
static Real LP2Di_(const Real *xarr, const Real *yarr,
                   const int lenx, const int ix, const int iy,
                   const Real *data, const Real x, const Real y);


/*----------------------------------------------------------------------------*/
/* IMPLEMENTATION of FUCNTIONS                                                */
/*----------------------------------------------------------------------------*/

Real fH2(const Real nH, const Real T, const Real Z_d, const Real xi_CR, 
         const Real G_H2) {
  Real kgr = 3.0e-17*Z_d;
  Real R = kgr * nH / (2. * xi_CR);
  Real a = 1.65*0.7;
  Real c = R;
  Real k_FUV = 5.7e-11 * G_H2;
  Real b = - (1.65*1.5 + 2.*R + k_FUV/(2.*xi_CR));
  Real x_H2 = (-b - sqrt(b*b - 4*a*c) )/(2.*a);
  return x_H2;
}

Real fCplus(const Real x_e, const Real x_H2, const Real nH, const Real T,
            const Real Z_d, const Real Z_g, const Real xi_CR, 
            const Real G_PE, const Real G_CI) {
  const Real xCtot=xCstd*Z_g;
  const Real small_ = 1e-50;
  Real k_C_cr = 3.85 * xi_CR;
  Real k_C_photo = 3.5e-10*G_CI;
  Real k_Cplus_e = CII_rec_rate_(T);
  Real psi_gr = 1.7 * G_PE * sqrt(T)/(nH * x_e + small_) + small_;
  const Real cCp_[7] = {45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2,
                        0.8120, 1.333e-4};
  Real k_Cplus_gr = 1.0e-14 * cCp_[0] / 
		           (
			           1.0 + cCp_[1]*pow(psi_gr, cCp_[2]) * 
								   (1.0 + cCp_[3] * pow(T, cCp_[4])
										             *pow( psi_gr, -cCp_[5]-cCp_[6]*log(T) ) 
									 ) 
								) * Z_d;
  Real k_Cplus_H2 = 3.3e-13 * pow(T, -1.3) * exp(-23./T);
  Real c = (k_C_cr + k_C_photo) / nH;
  Real al = k_Cplus_e*x_e + k_Cplus_gr + k_Cplus_H2*x_H2 + c;
  Real ar = xCtot * c;
  Real x_Cplus = ar / al;
  return x_Cplus;
}

Real fHplus(const Real x_e, const Real x_Cplus, const Real x_H2,
            const Real nH, const Real T, const Real Z_d, const Real xi_CR, 
            const Real G_PE) {
  const Real small_ = 1e-50;
  Real x_H = 1. - 2. * x_H2 - (x_e - x_Cplus);
  x_H = MAX(x_H, 0.0);
  Real k_Hplus_e = 2.753e-14 * pow( 315614.0 / T, 1.5) * pow( 
               1.0 + pow( 115188.0 / T, 0.407) , -2.242 );
  const Real cHp_[7] = {12.25, 8.074e-6, 1.378, 5.087e2,
                               1.586e-2, 0.4723, 1.102e-5};
  Real psi_gr = 1.7 * G_PE * sqrt(T)/(nH * x_e + small_) + small_;
  Real k_Hplus_gr = 1.0e-14 * cHp_[0] / 
             (
               1.0 + cHp_[1]*pow(psi_gr, cHp_[2]) * 
                 (1.0 + cHp_[3] * pow(T, cHp_[4])
                               *pow( psi_gr, -cHp_[5]-cHp_[6]*log(T) ) 
                 ) 
              ) * Z_d;
  Real k_cr_H = xi_CR * (2.3*x_H2 + 1.5*x_H);
  Real k_coll = 0;
  Real lnTe = log(T * 8.6173e-5);
  const Real T_coll = 7.0e2;
  if (T > T_coll) {
    k_coll = exp( -3.271396786e1 + 
                        (1.35365560e1 + (- 5.73932875 + (1.56315498 
                      + (- 2.877056e-1 + (3.48255977e-2 + (- 2.63197617e-3
                      + (1.11954395e-4 + (-2.03914985e-6)
          *lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe)*lnTe
                       );
  }
  Real c = k_cr_H * x_H + k_coll*x_e*nH;
  Real x_Hplus = c/( nH*(k_Hplus_e *  x_e + k_Hplus_gr) );
  x_Hplus = MIN(x_Hplus, 1.0);
  return x_Hplus;
}

Real fions(const Real x_e, const Real x_H2, const Real nH, const Real T,
           const Real Z_d, const Real Z_g, const Real xi_CR, 
           const Real G_PE, const Real G_CI) {
  const Real x_Cplus = fCplus(x_e, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI);
  const Real x_Hplus = fHplus(x_e, x_Cplus, x_H2, nH, T, Z_d, xi_CR, G_PE);
  return (x_Cplus + x_Hplus);
}

Real fe(const Real x_H2, const Real nH, const Real T,
        const Real Z_d, const Real Z_g, const Real xi_CR, 
        const Real G_PE, const Real G_CI) { 
  const Real rtol = 5e-2;
  const Real small_x = 1e-6;
  const Real small_f = 1e-20;
  const int maxiter = 20;
  Real x = 0;
  Real f = 0;
  Real xnew = 0;
  Real fnew = 0;
  int niter = 0;
  Real xprev = 0.5;
  Real fprev = 0;
  Real a = 0.0;
  Real fa = 0.0;
  Real b = 1.0;
  bool flag = true;
  niter = 0;
  fprev = fions(xprev, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI) - xprev;
  while (1) {
    x = fions(xprev, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI);
    if ( fabs((x-xprev)/(xprev+small_x)) < rtol || fabs(fprev) < small_x || 
         fabs((a-b)/(b+small_x)) < rtol) {
      break;
    }
    if (niter > maxiter) {
      printf("fe(): WARNING: niter>maxiter(=%d), x=%.2e, xprev=%.2e\n", 
             maxiter, x, xprev);
      printf("a=%.2e, b=%.2e\n", a, b);
      //TODO: throw ath_error
      break;
    }
    f = fions(x, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI) - x;
    if ( fabs(f - fprev) < small_f ) {
      xnew = (x + xprev)/2.;
    } else {
      xnew = x - f * (x-xprev)/(f-fprev);
    }
    if (xnew < a || xnew > b) {
      xnew = (a + b)/2.;
    }
    if (flag) {
      fa = fions(a, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)-a;
    }
    fnew = fions(xnew, x_H2, nH, T, Z_d, Z_g, xi_CR, G_PE, G_CI)-xnew;
    if (fa * fnew < 0.0) {
      b = xnew;
      flag = false;
    } else {
      a = xnew;
      flag = true;
    }
    xprev = xnew;
    fprev = fnew;
    niter = niter + 1;
  }
  return xprev;
}

Real fCO(const Real x_H2, const Real x_Cplus, const Real nH, 
         const Real Z_d, const Real Z_g, const Real xi_CR, const Real G_CO) {
  Real xCtot = xCstd * Z_g;
  Real kcr16 = xi_CR/1.0e-16;
  Real term1 = 4.0e3*Z_d/(kcr16*kcr16);
  Real ncrit2 = 0;
  Real x_CO = 0;
  Real x_CO_max1 = 0;
  Real x_CO_max2 = 0;
  term1 = MAX(term1, 1.);
  ncrit2 = 2.*pow(term1, pow(G_CO, 1./3.)) * (50.*kcr16/pow(Z_d, 1.4));
  if (nH >= ncrit2) {
    x_CO = 1.;
  } else {
    x_CO = nH/ncrit2;
  }
  x_CO *= xCtot;
  x_CO_max1 = xCtot - x_Cplus;
  x_CO_max2 = xCtot * x_H2 * 2.;
  x_CO = MIN(x_CO, x_CO_max1);
  x_CO = MIN(x_CO, x_CO_max2);
  return x_CO;
}

Real heatingCR(const Real x_e, const Real x_HI, const Real x_H2,
              const Real nH, const Real xi_CR) {
	/* ionization rate per H*/
  const Real kcr_H_fac = 1.15 * 2*x_H2 + 1.5 * x_HI;
  const Real kHI = xi_CR * kcr_H_fac;
  const Real kHe = xi_CR * 1.1;
  const Real kH2 = xi_CR * 2.0 * kcr_H_fac;
	const Real ktot = kHI*x_HI + kHe*xHetot + kH2*x_H2;
	/* heating rate per ionization in atomic region. 
	 * Draine ISM book eq (30.1)*/
  Real qHI;
  if (x_e > 1.0e-9) {
	  qHI = ( 6.5 + 26.4 * sqrt( x_e / (x_e+0.07) ) ) * eV_;
  } else { //prevent sqrt of small negative number
	  qHI =  6.5 * eV_;
  }

	/* Heating rate per ioniztion in molecular region.
	 * Despotic paper Appendix B*/
	Real qH2;
  const Real lognH = log10(nH);
  if (nH < 100.) { //prevent log of small negative number
    qH2 = 10. * eV_;
  } else if (lognH < 4) {
    qH2 = ( 10. + 3.*(lognH - 2.)/2. ) * eV_;
  } else if (lognH < 7) {
    qH2 = ( 13. + 4.*(lognH - 4.)/3. ) * eV_;
  } else if (lognH < 10) {
    qH2 = ( 17. + (lognH - 7.)/3. ) * eV_;
  } else {
    qH2 = 18. * eV_;
  }
	const Real qtot = x_HI*qHI + 2*x_H2*qH2;
	return (ktot*qtot);
}

Real heatingPE(const Real x_e, const Real nH, const Real T, 
               const Real Z_d, const Real G_PE) {
  const Real ne = x_e * nH;
  const Real x = 1.7 * G_PE * sqrt(T)/ne + 50.;
  const Real fac = ( CPE_[0] + CPE_[1]*pow(T, CPE_[4]) ) / 
    (
     1. + CPE_[2]*pow(x, CPE_[5]) * ( 1. + CPE_[3]*pow(x, CPE_[6]) ) 
     );
  const Real heating = 1.7e-26 * G_PE * Z_d * fac;
  return heating;
}

Real heatingH2pump(const Real x_HI, const Real x_H2, const Real nH,
                   const Real T, const Real G_H2) {
  const Real dot_xH2_photo = 5.7e-11 * G_H2 * x_H2;
  const Real de = 1.6 * x_HI * exp( - pow(400./T, 2) )
                    + 1.4 * x_H2 * exp( - (12000./(T + 1200.)));
  const Real ncr = 1.0e6 / sqrt(T) / de;
  const Real f = 1. / (1. + ncr/nH);
  return dot_xH2_photo * 9. * 2.2*f * eV_;
}

Real CII_rec_rate_(const Real temp) {
  Real A, B, T0, T1, C, T2, BN, term1, term2, alpharr, alphadr;
  A = 2.995e-9;
  B = 0.7849;
  T0 =  6.670e-3;
  T1 = 1.943e6;
  C = 0.1597;
  T2 = 4.955e4;
  BN = B + C * exp(-T2/temp);
  term1 = sqrt(temp/T0);
  term2 = sqrt(temp/T1);
  alpharr = A / ( term1*pow(1.0+term1, 1.0-BN) * pow(1.0+term2, 1.0+BN) );
  alphadr = pow( temp, -3.0/2.0 ) * ( 6.346e-9 * exp(-1.217e1/temp) +
        9.793e-09 * exp(-7.38e1/temp) + 1.634e-06 * exp(-1.523e+04/temp) );
  return (alpharr+alphadr);
}

Real q10CII_(const Real nHI, const Real nH2, const Real ne, const Real T) {
	/*Draine (2011) ISM book eq (17.16) and (17.17)*/
	const Real T2 = T/100.;
	const Real k10e = 4.53e-8 * sqrt(1.0e4/T);
	const Real k10HI = 7.58e-10 * pow(T2, 0.1281+0.0087*log(T2));
	Real k10oH2 = 0;
	Real k10pH2 = 0;
	Real tmp = 0;
	if (T < 500.) {
    /*fit in Wiesenfeld & Goldsmith 2014*/
    k10oH2 = (5.33 + 0.11*T2)*1.0e-10;
    k10pH2 = (4.43 + 0.33*T2)*1.0e-10;
  } else {
    /* Glover+ Jappsen 2007, for high temperature scales similar to HI*/
		tmp = pow(T, 0.07);
		k10oH2 = 3.74757785025e-10*tmp;
		k10pH2 = 3.88997286356e-10*tmp;
	}
	const Real k10H2 = k10oH2*fo_ + k10pH2*fp_;
	//printf("q10e=%0.4e, q10HI=%0.4e, q10H2=%0.4e\n", k10e*ne, k10HI*nHI, k10H2*nH2);
	return (k10e*ne + k10HI*nHI + k10H2*nH2);
}

Real cooling2Level_(const Real q01, const Real q10, const Real A10,
                    const Real E10, const Real xs) {
	const Real f1 = q01 / (q01 + q10 + A10);
	return f1*A10*E10*xs;
}

Real cooling3Level_(const Real q01, const Real q10, const Real q02,
                    const Real q20, const Real q12, const Real q21,
										const Real A10, const Real A20, const Real A21,
                    const Real E10, const Real E20, const Real E21,
										const Real xs) {
	const Real R10 = q10 + A10;
	const Real R20 = q20 + A20;
	const Real R21 = q21 + A21;
	const Real a0 = R10*R20 + R10*R21 + q12*R20;
	const Real a1 = q01*R20 + q01*R21 + R21*q02;
	const Real a2 = q02*R10 + q02*q12 + q12*q01;
	const Real de = a0 + a1 + a2;
	const Real f1 = a1 / de;
	const Real f2 = a2 / de;
	return ( f1*A10*E10 + f2*(A20*E20 + A21*E21) )*xs;
}

int linearInterpIndex_(const int len, const Real xarr[], const Real x){
  int i = 0;
  if ( x < xarr[0]) {
    return 0;
  } else if ( x > xarr[len-1]) {
    return len-2;
  } else {
    for (i=0; x>xarr[i]; i++) {}
    return i-1;
  }
}

Real linearInterp_(const Real x0, const Real x1, const Real y0,
                   const Real y1, const Real x){
  return y0 + ( (y1-y0)/(x1-x0) ) * (x-x0);
}

Real LP1Di_(const Real *xarr, const Real *data, const int ix,
            const Real x) {
  return linearInterp_(xarr[ix], xarr[ix+1], data[ix], data[ix+1], x);
}

Real LP2Di_(const Real *xarr, const Real *yarr,
            const int lenx, const int ix, const int iy,
            const Real *data, const Real x, const Real y) {
  Real fl1, fl2;
  const Real x0 = xarr[ix];
  const Real x1 = xarr[ix+1];
  fl1 = linearInterp_(x0, x1, data[iy*lenx + ix], data[iy*lenx + ix+1], x);
  fl2 = linearInterp_(x0, x1, data[(iy+1)*lenx + ix], data[(iy+1)*lenx + ix+1], x);
  return linearInterp_(yarr[iy], yarr[iy+1], fl1, fl2, y);
}

Real coolingLya(const Real x_e, const Real x_HI, const Real nH, const Real T) {
  const Real ne = x_e * nH;
  const Real A = 6.3803e-9;
  const Real beta = 1.17;
  const Real T4 = T / 1.0e4;
  const Real fac = A * pow(T4, beta);
  const Real k01e = fac * exp(-11.84/T4);
  const Real q01 = k01e * ne;
  const Real q10 = (g0HI_/g1HI_) * fac * ne;
	return cooling2Level_(q01, q10, A10HI_, E10HI_, x_HI);
}

Real coolingOI(const Real x_e, const Real x_OI, const Real x_HI,
               const Real x_H2, const Real nH, const Real T) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
  /*collisional rates from  Draine (2011) ISM book Appendix F Table F.6*/
  const Real T2 = T/100;
  const Real lnT2 = log(T2);
  /*HI*/
  const Real k10HI = 3.57e-10 * pow(T2, 0.419-0.003*lnT2);
  const Real k20HI = 3.19e-10 * pow(T2, 0.369-0.006*lnT2);
  const Real k21HI = 4.34e-10 * pow(T2, 0.755-0.160*lnT2);
  /*H2*/
  const Real k10H2p = 1.49e-10 * pow(T2, 0.264+0.025*lnT2);
  const Real k10H2o = 1.37e-10 * pow(T2, 0.296+0.043*lnT2);
  const Real k20H2p = 1.90e-10 * pow(T2, 0.203+0.041*lnT2);
  const Real k20H2o = 2.23e-10 * pow(T2, 0.237+0.058*lnT2);
  const Real k21H2p = 2.10e-12 * pow(T2, 0.889+0.043*lnT2);
  const Real k21H2o = 3.00e-12 * pow(T2, 1.198+0.525*lnT2);
  const Real k10H2 = k10H2p*fp_ + k10H2o*fo_;
	const Real k20H2 = k20H2p*fp_ + k20H2o*fo_;
	const Real k21H2 = k21H2p*fp_ + k21H2o*fo_;
  /*e*/
  /*fit from Bell+1998*/
  const Real k10e = 5.12e-10 * pow(T, -0.075);
  const Real k20e = 4.86e-10 * pow(T, -0.026);
  const Real k21e = 1.08e-14 * pow(T, 0.926);
  /*total collisional rates*/
	const Real q10 = k10HI*nHI + k10H2*nH2 + k10e * ne;
	const Real q20 = k20HI*nHI + k20H2*nH2 + k20e * ne;
	const Real q21 = k21HI*nHI + k21H2*nH2 + k21e * ne;
	const Real q01 = (g1OI_/g0OI_) * q10 * exp( -E10OI_/(kb_*T) );
	const Real q02 = (g2OI_/g0OI_) * q20 * exp( -E20OI_/(kb_*T) );
	const Real q12 = (g2OI_/g1OI_) * q21 * exp( -E21OI_/(kb_*T) );

	return cooling3Level_(q01,q10, q02, q20, q12, q21, A10OI_,A20OI_,
												A21OI_, E10OI_, E20OI_, E21OI_, x_OI);
}

Real coolingCII(const Real x_e, const Real x_Cplus, const Real x_HI,
                const Real x_H2, const Real nH, const Real T) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
	const Real q10 = q10CII_(nHI, nH2, ne, T);
	const Real q01 = (g1CII_/g0CII_) * q10 * exp( -E10CII_/(kb_*T) );
	return cooling2Level_(q01, q10, A10CII_, E10CII_, x_Cplus);
}

Real coolingCI(const Real x_e, const Real x_CI, const Real x_HI,
               const Real x_H2, const Real nH, const Real T) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
	/*e collisional coefficents from Johnson, Burke, & Kingston 1987, 
	 * JPhysB, 20, 2553*/
	const Real T2 = T/100.;
	const Real lnT2 = log(T2);
	const Real lnT = log(T);
	/*ke(u,l) = fac*gamma(u,l)/g(u)*/
	const Real fac = 8.629e-8 * sqrt(1.0e4/T);
	Real k10e, k20e, k21e;
	Real lngamma10e, lngamma20e, lngamma21e; /*collisional strength*/
	if (T < 1.0e3) {
		lngamma10e = (((-6.56325e-4*lnT -1.50892e-2)*lnT + 3.61184e-1)*lnT 
				          -7.73782e-1)*lnT - 9.25141;
		lngamma20e = (((0.705277e-2*lnT - 0.111338)*lnT +0.697638)*lnT
				          - 1.30743)*lnT -7.69735;
		lngamma21e = (((2.35272e-3*lnT - 4.18166e-2)*lnT +0.358264)*lnT 
				          - 0.57443)*lnT -7.4387;

	} else {
		lngamma10e = (((1.0508e-1*lnT - 3.47620)*lnT + 4.2595e1)*lnT
									- 2.27913e2)*lnT + 4.446e2;
		lngamma20e = (((9.38138e-2*lnT - 3.03283)*lnT +3.61803e1)*lnT 
									- 1.87474e2)*lnT +3.50609e2;
		lngamma21e = (((9.78573e-2*lnT - 3.19268)*lnT +3.85049e1)*lnT 
				          - 2.02193e2)*lnT +3.86186e2;
	}
	k10e = fac * exp(lngamma10e) / g1CI_;
	k20e = fac * exp(lngamma20e) / g2CI_;
	k21e = fac * exp(lngamma21e) / g2CI_;
	/*HI collisional rates, Draine (2011) ISM book Appendix F Table F.6
	 * NOTE: this is more updated than the LAMBDA database.*/
	const Real k10HI = 1.26e-10 * pow(T2, 0.115+0.057*lnT2);
	const Real k20HI = 0.89e-10 * pow(T2, 0.228+0.046*lnT2);
	const Real k21HI = 2.64e-10 * pow(T2, 0.231+0.046*lnT2);
	/*H2 collisional rates, Draine (2011) ISM book Appendix F Table F.6*/
	const Real k10H2p = 0.67e-10 * pow(T2, -0.085+0.102*lnT2);
	const Real k10H2o = 0.71e-10 * pow(T2, -0.004+0.049*lnT2);
	const Real k20H2p = 0.86e-10 * pow(T2, -0.010+0.048*lnT2);
	const Real k20H2o = 0.69e-10 * pow(T2, 0.169+0.038*lnT2);
	const Real k21H2p = 1.75e-10 * pow(T2, 0.072+0.064*lnT2);
	const Real k21H2o = 1.48e-10 * pow(T2, 0.263+0.031*lnT2);
	const Real k10H2 = k10H2p*fp_ + k10H2o*fo_;
	const Real k20H2 = k20H2p*fp_ + k20H2o*fo_;
	const Real k21H2 = k21H2p*fp_ + k21H2o*fo_;
	/* The totol collisonal rates*/
	const Real q10 = k10HI*nHI + k10H2*nH2 + k10e*ne;
	const Real q20 = k20HI*nHI + k20H2*nH2 + k20e*ne;
	const Real q21 = k21HI*nHI + k21H2*nH2 + k21e*ne;
	const Real q01 = (g1CI_/g0CI_) * q10 * exp( -E10CI_/(kb_*T) );
	const Real q02 = (g2CI_/g0CI_) * q20 * exp( -E20CI_/(kb_*T) );
	const Real q12 = (g2CI_/g1CI_) * q21 * exp( -E21CI_/(kb_*T) );
	return cooling3Level_(q01,q10, q02, q20, q12, q21, A10CI_,A20CI_,
												A21CI_, E10CI_, E20CI_, E21CI_, x_CI);
}

Real coolingCO(const Real x_e, const Real x_CO, const Real x_HI,
               const Real x_H2, const Real nH, const Real T, 
               const Real dvdr) {
  const Real nHI = x_HI * nH;
  const Real nH2 = x_H2 * nH;
  const Real ne = x_e * nH;
  const Real nCO = x_CO * nH;
  //effective column of CO
  //maximum escape probability length, in cgs unites
  const Real Leff_CO_max = 3.086e20; //100 pc 
  const Real mCO = 4.68e-23;
  const Real vth = sqrt(2. * kb_ * T / mCO);
	const Real grad_small = vth/Leff_CO_max;
  const Real gradeff = MAX(dvdr, grad_small);
  const Real NCOeff = nCO / gradeff;
  //maximum temperature above which use Tmax for cooling rate interpolation
  const Real Tmax_CO = 2000.; 
  Real T1 = 0;;
  if (T < Tmax_CO) {
    T1 = T;
  } else {
    T1 = Tmax_CO;
  }
  const Real facT = pow(1. - exp(-T1), 1.0e3);
  /*small number for a very small NCOeff*/
  const Real eps = 1.0e13;
  const Real log_NCOeff = log10(NCOeff*1.0e5 + eps); /*unit: cm^-2 / (km/s) */
  const Real Troot4 = pow(T1, 0.25);
  const Real neff = nH2 + 1.75*Troot4 * nHI + 680.1/Troot4 * ne;
  /* interpolate parameters using given T and NCOeff*/
  /* index of T and Neff*/
  const int iT0 =  linearInterpIndex_(lenTCO_, TCO_, T1);
  const int iNeff0 = linearInterpIndex_(lenNeffCO_, NeffCO_, log_NCOeff);
  /* L0 */
  const Real log_L0 = - LP1Di_(TCO_, L0CO_, iT0, T1);
  const Real L0 = pow(10, log_L0);
  /* LLTE */
  const Real log_LLTE = - LP2Di_(TCO_, NeffCO_, lenTCO_, iT0, iNeff0, 
                                           LLTECO_, T1, log_NCOeff);
  const Real LLTE = pow(10, log_LLTE);
  /* n1/2*/
  const Real log_nhalf = LP2Di_(TCO_, NeffCO_, lenTCO_, iT0, iNeff0, 
                                          nhalfCO_, T1, log_NCOeff);
  const Real nhalf = pow(10, log_nhalf);
  /* alpha*/
  const Real alpha = LP2Di_(TCO_, NeffCO_, lenTCO_, iT0, iNeff0, 
                                      alphaCO_, T1, log_NCOeff);
  const Real inv_LCO = 1./L0 + neff/LLTE 
                         + 1./L0 * pow(neff/nhalf, alpha) * (1. - nhalf*L0/LLTE);
  return (1./inv_LCO) * neff * x_CO * facT;
}

Real coolingRec(const Real x_e, const Real nH, const Real T, 
                const Real Z_d, const Real G_PE) {
  const Real ne = x_e * nH;
  const Real x = 1.7 * G_PE * sqrt(T)/ne + 50.;
  const Real lnx = log(x);
  const Real cooling = 1.0e-28 * ne * pow(T, DPE_[0] + DPE_[1]/lnx) 
                          * exp( DPE_[2] + (DPE_[3] - DPE_[4]*lnx)*lnx );
  return cooling * Z_d;
}

Real coolingHot(const Real T, const Real Z_g) {
  const Real log10T = log10(T);
  const int iT =  linearInterpIndex_(len_hot_, log10T_hot_, log10T);
  const Real log10L_HHe = LP1Di_(log10T_hot_, log10L_hot_HHe_, iT, log10T);
  const Real log10L_metal = LP1Di_(log10T_hot_, log10L_hot_metal_, iT, log10T);
  const Real L_HHe = pow(10, log10L_HHe);
  const Real L_metal = pow(10, log10L_metal) * Z_g;
  return (L_HHe + L_metal);
}

void get_abundances(const Real nH, const Real T, const Real dvdr, const Real Z,
                    const Real xi_CR, const Real G_PE, const Real G_CI,
                    const Real G_CO, const Real G_H2,
                    Real *px_e, Real *px_HI, Real *px_H2, Real *px_Cplus,
                    Real *px_CI, Real *px_CO, Real *px_OI) {
  Real x_e, x_HI, x_H2, x_Hplus, x_Cplus, x_CI, x_CO, x_OI;
  x_H2 = fH2(nH, T, Z, xi_CR, G_H2);
  x_e = fe(x_H2, nH, T, Z, Z, xi_CR, G_PE, G_CI);
  x_Cplus = fCplus(x_e, x_H2, nH, T, Z, Z, xi_CR, G_PE, G_CI);
  x_Hplus = fHplus(x_e, x_H2, x_Cplus, nH, T, Z, xi_CR, G_PE);
  x_CO = fCO(x_H2, x_Cplus, nH, Z, Z, xi_CR, G_CO);
  x_HI = MAX(1. - 2.0*x_H2 - x_Hplus, 0.);
  x_CI = MAX(xCstd*Z - x_CO - x_Cplus, 0.);
  x_OI = MAX(xOstd*Z - x_CO, 0.);
  *px_e = x_e;
  *px_HI = x_HI;
  *px_H2 = x_H2;
  *px_Cplus=x_Cplus;
  *px_CI = x_CI;
  *px_CO = x_CO;
  *px_OI = x_OI;
  return;
}

Real heating(const Real x_e, const Real x_HI, const Real x_H2,
             const Real nH, const Real T, const Real Z,
             const Real xi_CR, const Real G_PE, const Real G_H2) {
  const Real h_PE = heatingPE(x_e, nH, T, Z, G_PE);
  const Real h_CR = heatingCR(x_e, x_HI, x_H2, nH, xi_CR);
  const Real h_H2pump = heatingH2pump(x_HI, x_H2, nH, T, G_H2);
  return (h_PE + h_CR + h_H2pump)/nH;
}

Real cooling(const Real x_e, const Real x_HI, const Real x_H2,
             const Real x_Cplus, const Real x_CI, const Real x_CO, const Real x_OI, 
             const Real nH, const Real T, const Real dvdr,
             const Real Z, const Real G_PE) {
  const Real T_hot = 1.2e4;
  Real c_Lya, c_OI, c_CII, c_CI, c_CO, c_Rec;
  if (T > T_hot) {
    return coolingHot(T, Z);
  } else {
    c_Lya = coolingLya(x_e, x_HI, nH, T);
    c_OI = coolingOI(x_e, x_OI, x_HI, x_H2, nH, T);
    c_CII = coolingCII(x_e, x_Cplus, x_HI, x_H2, nH, T);
    c_CI = coolingCI(x_e, x_CI, x_HI, x_H2, nH, T);
    c_CO = coolingCO(x_e, x_CO, x_HI, x_H2, nH, T, dvdr);
    c_Rec = coolingRec(x_e, nH, T, Z, G_PE);
    return (c_Lya + c_OI + c_CII + c_CI + c_CO + c_Rec)/nH;
  }
}
