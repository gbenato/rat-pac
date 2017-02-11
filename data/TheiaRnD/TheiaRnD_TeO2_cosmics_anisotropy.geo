/////////////////////
// DARK BOX
/////////////////////

{
  name: "GEO",
  index: "world",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 1, // omitted for visualization
  mother: "",
  type: "box",
  size: [3000.0,3000.0,600.0], //mm, half-lenght
  material: "air",
}

{
  name: "GEO",
  index: "darkbox",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 1, // omitted for visualization
  mother: "world",
  type: "box",
  size: [762.0,762.0,508.0], //mm, half-lenght
  material: "acrylic_black", //acrylic_black
  color: [0.5, 0.2, 0.1, 0.1],
}

{
  name: "GEO",
  index: "inner",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 1, // omitted for visualization
  mother: "darkbox",
  type: "box",
  size: [711.2,711.2,457.2], //mm, half-lenght
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}

/////////////////////
// ACRYLIC BLOCK
/////////////////////
{
  name: "GEO",
  index: "block",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  size: [160.0,160.0,32.5],
  position: [-398.0, -367.0, -270.0],
//  size: [160.0,160.0,77.5], //alternative box as mothervolume around the entire TeO2 (technically avoids pyramids...)
//  position: [-398.0, -367.0, -235.0], // 45 mm higher center, 90 mm thicker can enclose Teo2 standing on edge 
  material: "chsrc_uvt_acrylic",
//  material: "tellurium_dioxyde",
//  material: "air",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "block1",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  size: [75.0,25.0,25.0],
  position: [-498.0, -367.0, -212.5],
  material: "air",
  color: [0.1, 0.3, 0.8, 0.1],
}


{
  name: "GEO",
  index: "block2",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  size: [75.0,25.0,25.0],
  position: [-298.0, -367.0, -212.5],
  material: "air",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "block3",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  size: [100.0,60,25.0],
  position: [-398.0, -452.0, -212.5],
  material: "air",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "block4",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  size: [100.0,60.0,25.0],
  position: [-398.0, -282.0, -212.5],
  material: "air",
  color: [0.1, 0.3, 0.8, 0.1],
}

/////////////////////
// COVER
/////////////////////
{
  name: "GEO",
  index: "cover",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "block",
  type: "box",
  size: [160.0,160.0,0.5],
  position: [0.0,0.0,32.0],
  material: "chsrc_uvt_acrylic",
  color: [0.5, 0.2, 0.1, 0.1],
//  surface: "mirror",
}

{
  name: "GEO",
  index: "cover_window",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "cover",
  type: "tube",
  r_max: 50.0,
  size_z: 0.5,
  position: [0.0,0.0,0.0],
  material: "chsrc_uvt_acrylic",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "cover_hollow",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "cover_window",
  type: "tube",
  r_max: 5.0,
  size_z: 0.5,
  position: [0.0,0.0,0.0],
  material: "acrylic_black",
  color: [0.0, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "hollow",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "block",
  type: "tube",
  position: [0.0, 0.0, -0.5],
  r_max: 5.0,
  size_z: 32.0,
  material: "acrylic_black", // acrylic_black
  color: [0.0, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "chip",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "block",
  type: "tube",
  position: [0.0, 0.0, -31.5],
  r_max: 10.0,
  r_min: 5.0,
  size_z: 1.0,
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}

/////////////////
//  TeO2 CRYSTAL, and indivdual surface volumes
/////////////////

{
  name: "GEO",
  index: "TeO2",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  size: [25.0,25.0,25.0], //mm, half-lenght
  position: [-398.0, -367.0, -212.5],
  rotation: [0.0, 0.0, 0.0]
//  material: "chsrc_uvt_acrylic"
  material: "tellurium_dioxyde", //tellurium_dioxyde
//  material: "air", //tellurium_dioxyde
  color: [0.8, 0.8, 0.8, 0.1],
}

{
  name: "GEO",
  index: "TeO2_bottom",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "TeO2",
  type: "box",
  size: [25.0,25.0,0.1], //mm, half-lenght
  position: [0.0, 0.0, -24.9],
  rotation: [0.0, 0.0, 0.0]
//  material: "chsrc_uvt_acrylic"
  material: "tellurium_dioxyde", //tellurium_dioxyde
  color: [0.8, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "TeO2_top", // surface volume that allos to specify specific surface border properties
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "TeO2",
  type: "box",
  size: [25.0,25.0,0.1], //mm, half-lenght
  position: [0.0, 0.0, 24.9],
  rotation: [0.0, 0.0, 0.0]
//  material: "chsrc_uvt_acrylic"
//  material: "tellurium_dioxyde", //tellurium_dioxyde
  material: "air", //tellurium_dioxyde
  color: [0.8, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "TeO2_posX",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "TeO2",
  type: "box",
  size: [0.1,25.0,24.8], //mm, half-lenght
  position: [24.9, 0.0, 0.0],
  rotation: [0.0, 0.0, 0.0]
//  material: "chsrc_uvt_acrylic"
//  material: "tellurium_dioxyde", //tellurium_dioxyde
  material: "air", //tellurium_dioxyde
  color: [0.8, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "TeO2_negX",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "TeO2",
  type: "box",
  size: [0.1,25.0,24.8], //mm, half-lenght
  position: [-24.9, 0.0, 0.0],
  rotation: [0.0, 0.0, 0.0]
//  material: "chsrc_uvt_acrylic"
//  material: "tellurium_dioxyde", //tellurium_dioxyde
  material: "air", //tellurium_dioxyde
  color: [0.8, 0.0, 0.0, 0.1],
}

{  
  name: "GEO",
  index: "TeO2_posY",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "TeO2",
  type: "box",
  size: [24.8,.1,24.8], //mm, half-lenght
  position: [0.0, 24.9, 0.0],
  rotation: [0.0, 0.0, 0.0]
//  material: "chsrc_uvt_acrylic"
//  material: "tellurium_dioxyde", //tellurium_dioxyde
  material: "air", //tellurium_dioxyde
  color: [0.8, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "TeO2_negY",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "TeO2",
  type: "box",
  size: [24.8,0.1,24.8], //mm, half-lenght
  position: [0.0, 24.9, 0.0],
  rotation: [0.0, 0.0, 0.0]
//  material: "chsrc_uvt_acrylic"
//  material: "tellurium_dioxyde", //tellurium_dioxyde
  material: "air", //tellurium_dioxyde
  color: [0.8, 0.0, 0.0, 0.1],
}

/////////////////
//  VESSEL
/////////////////
//{
//  name: "GEO",
//  index: "vessel",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "inner",
//  type: "CheSSVessel",
//  r_max: 50.0,
//  size_z: 15.0,
//  rotation:  [0.0, 0.0, 135.0],
//  position: [-398.0, -367.0, -220.91],
//  material: "chsrc_uvt_acrylic",
//  color: [0.1, 0.3, 0.8, 0.1],
//}

//{
//  name: "GEO",
//  index: "cavity",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "vessel",
//  type: "tube",
//  r_max: 20.0,
//  size_z: 1.905,
//  position: [0.0, 0.0, 14.685],
//  material: "air",
//  color: [0.0, 0.0, 0.0, 0.1],
//}

//{
//  name: "GEO",
//  index: "content",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "vessel",
//  type: "tube",
//  r_max: 45.0, //LAB->45mm, WBLS->35mm
//  size_z: 12.89,
//  position: [0.0, 0.0, -1.7],
//  material: "water",
//  color: [0.5, 0.1, 0.5, 0.5],
//}

/////////////////////////////
// RADIACTIVE SOURCES
/////////////////////////////
//{
//  name: "GEO",
//  index: "envelope",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "inner",
//  type: "tube",
//  position: [-398.0, -367.0, -185.9],
//  rotation:  [0.0, 0.0, 0.0],
//  r_max: 12.808,
//  size_z: 1.587, //half height
//  material: "acrylic_black",
//  color: [0.1, 1.0, 0.3, 0.8],
//}

//{
//  name: "GEO",
//  index: "source",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "envelope",
//  type: "tube",
//  position: [0.0, 0.0, -0.254],
//  r_max: 3.175,
//  size_z: 1.333, //half height
//  material: "strontium", //strontium
//  color: [0.1, 1.0, 1.0, 0.8],
//}

//{
//  name: "GEO",
//  index: "envelope",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "inner",
//  type: "tube",
//  position: [-408.76, -417.98, -222.5],
//  rotation:  [90.0, -11.914, 0.0],
//  r_max: 12.808,
//  size_z: 1.587, //half height
//  material: "acrylic_black",
//  color: [0.1, 1.0, 0.3, 0.8],
//}

//{
//  name: "GEO",
//  index: "source",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "envelope",
//  type: "tube",
//  position: [0.0, 0.0, 0.254],
//  r_max: 3.175,
//  size_z: 1.333, //half height
//  material: "strontium", //strontium
//  color: [0.1, 1.0, 1.0, 0.8],
//}

///////////////////////
// COSMIC TAGS
///////////////////////
{
  name: "GEO",
  index: "tag1",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "tube",
  position: [-398.0, -367.0, -147.987], //3.5 cm higher than usual 
  r_max: 5.0,
  size_z: 25.0,
  material: "acrylic_black",
  color: [0.0, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "tag2",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "tube",
  position: [-398.0, -367.0, -327.5],
  r_max: 5.0,
  size_z: 25.0,
  material: "acrylic_black",
  color: [0.0, 0.0, 0.0, 0.1],
}
/////////////////////////////////


//////////////////
// PMTS
//////////////////

//Why is this container outcommented - check the original geo!
// Container for pmts
//{
//  name: "GEO",
//  index: "pmt_holder",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "inner",
//  type: "box",
//  size: [110.0,110.0,15.0],
//  position: [-400.0, -400.0, 0.0],
//  material: "acrylic_black",
//  color: [0.5, 0.0, 0.0, 0.5],
//}

{
  name: "GEO",
  index: "ring_pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
  type: "pmtarray",
  pmt_model: "h11934", //h11934, 1-inch cubic
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_CROSS_SIDE",
  orientation: "manual",
}

{
  name: "GEO",
  index: "light_pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
  type: "pmtarray",
  pmt_model: "r7081_hqe", // big ones, for light yield measurement
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_CLOSE",
  orientation: "point",
  orient_point: [-398.0, -367.0, -237.0],
}

{
   name: "GEO",
   index: "trigger_pmt",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   mother: "inner",
   type: "pmtarray",
   pmt_model: "h11934",// again small, 1-inch cubic
   pmt_detector_type: "idpmt",
   sensitive_detector: "/mydet/pmt/inner",
   efficiency_correction: 1.0,
   pos_table: "PMTINFO_TRIGGER",
   orientation: "manual",
//   orient_point: [-398.0, -367.0, -222.2],
}

{
  name: "GEO",
  index: "muon_pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "world",
  type: "pmtarray",
  pmt_model: "h11934",
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_MUON_TAGS",
  orientation: "manual",
  orient_point: [-400.0, -400.0, -200.0],
}

{
   name: "GEO",
   index: "panels_pmt",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   mother: "world",
   type: "pmtarray",
   pmt_model: "h11934", //h11934, r7081_hqe, r11780_hqe, fast_test
   pmt_detector_type: "idpmt",
   sensitive_detector: "/mydet/pmt/inner",
   efficiency_correction: 1.0,
   pos_table: "PMTINFO_PANELS",
   orientation: "manual",
   orient_point: [0.0, 0.0, 400.0],
}


//////////////////
// Borders - Surfaces
//////////////////


{
   name: "GEO",
   index:"teo2top",
   mother: "TeO2",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   type: "border",
   volume1: "TeO2",
   volume2: "TeO2_top",
   surface: "teo2_soft_matt",
//   surface: "teo2_hard_glossy",
   reverse: 1.0,
}

{
   name: "GEO",
   index:"teo2posx",
   mother: "TeO2",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   type: "border",
   volume1: "TeO2",
   volume2: "TeO2_posX",
   surface: "teo2_soft_matt",
//   surface: "teo2_hard_glossy",
   reverse: 1.0,
}

{
   name: "GEO",
   index:"teo2negx",
   mother: "TeO2",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   type: "border",
   volume1: "TeO2",
   volume2: "TeO2_negX",
   surface: "teo2_soft_matt",
//   surface: "teo2_hard_glossy",
   reverse: 1.0,
}

{
   name: "GEO",
   index:"teo2posy",
   mother: "TeO2",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   type: "border",
   volume1: "TeO2",
   volume2: "TeO2_posY",
   surface: "teo2_soft_matt",
//   surface: "teo2_hard_glossy",
   reverse: 1.0,
}


{
   name: "GEO",
   index:"teo2negy",
   mother: "TeO2",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   type: "border",
   volume1: "TeO2",
   volume2: "TeO2_negY",
   surface: "teo2_soft_matt",
//   surface: "teo2_hard_glossy",
   reverse: 1.0,
}

{
   name: "GEO",
   index:"teo2bottom",
   mother: "TeO2",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   type: "border",
   volume1: "TeO2_bottom",
   volume2: "block",
//   surface: "teo2_acrylic_smooth",
   surface: "teo2_acrylic_rough",
   reverse: 1.0,
}
