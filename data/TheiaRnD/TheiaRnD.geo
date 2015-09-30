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
  invisible: 0, // omitted for visualization
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
  invisible: 0, // omitted for visualization
  mother: "darkbox",
  type: "box",
  size: [711.2,711.2,457.2], //mm, half-lenght
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}

///////////////////////////////////////
//             VESSELS               //
///////////////////////////////////////

///////////////////////////////////////
//// ACRYLIC BLOCK
{
  name: "GEO",
  index: "outer_vessel",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  size: [160.0,160.0,32.5],
  position: [-400.0, -400.0, -200.0],
  material: "acrylic_berkeley",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "cover",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "outer_vessel",
  type: "box",
  size: [160.0,160.0,0.5],
  position: [0.0,0.0,32.0],
//  position: [-400.0, -400.0, -180.0],
  material: "mirror",
  color: [0.5, 0.2, 0.1, 0.1],
}

{
  name: "GEO",
  index: "mirror",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "world", //not used
  type: "border",
  volume1: "inner",
  volume2: "cover",
  reverse: 1,
  surface: "mirror",
}

{
  name: "GEO",
  index: "cover_window",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "cover",
  //Cylinder
  type: "tube",
  r_max: 50.0,
  size_z: 0.5,
  position: [0.0,0.0,0.0],
  material: "acrylic_berkeley",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "cover_hollow",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "cover_window",
  //Cylinder
  type: "tube",
  r_max: 5.5,
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
  mother: "outer_vessel",
  type: "tube",
  position: [0.0, 0.0, -0.5],
  r_max: 5.5,
  size_z: 32.0,
  material: "air", // acrylic_black
  color: [0.0, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "container",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner", //inner_vessel
  //Cylinder
  type: "tube",
  r_max: 50.0,
  size_z: 15.0,
  //Cuboid
//  type: "box",
//  rotation:  [0.0, 0.0, 0.0],
//  size: [50.0,50.0,15.0],
  position: [-400.0, -400.0, -152.5],
  material: "acrylic_berkeley",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "content",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "container",
  type: "tube",
  r_max: 35.0, //LAB->35.0, WBLS->30.0
  size_z: 14.0,
  position: [0.0, 0.0, 1.0],
  material: "wbls_10pct", //water, wbls_1pct_mod, wbls_1pct, wbls_5pct, wbls_10pct, labppo_scintillator
  color: [0.5, 0.1, 0.5, 0.5],
}


// {
//   name: "GEO",
//   index: "blind_spot",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "container",
//   type: "tube",
//   position: [0.0, 0.0, -17.5], //(-20.0)
//   r_max: 10.0, //5.0
//   size_z: 2.5, //5.0
//   material: "acrylic_black", //air
//   color: [0.0, 0.0, 1.0, 1.0],
// }

//////////////////////////////////////////////

///////////////////////////////////////
//// BOX
// {
//   name: "GEO",
//   index: "vessel",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "inner",
//   type: "box",
//   position: [0.0, 0.0, 200.0],
// //  size: [20.0,20.0,12.0], //box
// //  size: [20.0,20.0,4.0], //thin box
// //  size: [10.0,10.0,20.0], //cuboid
//   size: [200.0,200.0,70.5], //tank
//   material: "quartz", //quartz, acrylic_berkeley
//   color: [0.1, 0.3, 0.8, 0.1],
// }

// {
//   name: "GEO",
//   index: "content",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "vessel",
//   type: "box",
//   position: [0.0, 0.0, 0.0],
//   size: [18.0,18.0,3.0], //mm, half-lenght
// //  size: [9.0,9.0,19.0], //mm, half-lenght
//   material: "labppo_scintillator",
//   color: [0.1, 0.1, 1.0, 0.5],
// }
//////////////////////////////////////////////

///////////////////////////////////////////////
// //// Cylinder
// {
//   name: "GEO",
//   index: "vessel",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "inner",
//   type: "tube",
//   r_max: 10.0, //mm
//   size_z: 20.0, //mm
//   position: [0.0, 0.0, 200.0],
//   material: "acrylic_berkeley",
//   color: [0.1, 0.3, 0.8, 0.1],
// }

// {
//   name: "GEO",
//   index: "content",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "vessel",
//   type: "tube",
//   position: [0.0, 0.0, 0.0],
//   r_max: 8.0, //mm
//   size_z: 18.0, //mm
//   material: "wbls_5pct", //labppo_scintillator, water
//   color: [0.1, 0.1, 1.0, 0.5],
// }
////////////////////////////////////////////////

//////////////////////////////////////////////
//// SPHERE
// {
//   name: "GEO",
//   index: "vessel",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "inner",
//   type: "sphere",
//   r_max: 20.0, //mm
//   position: [0.0, 0.0, 200.0],
//   material: "acrylic_berkeley",
//   color: [0.1, 0.3, 0.8, 0.1],
// }

// {
//   name: "GEO",
//   index: "content",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "vessel",
//   type: "sphere",
//   position: [0.0, 0.0, 0.0],
//   r_max: 18.0, //mm
//   material: "labppo_scintillator",
//   color: [0.1, 0.1, 1.0, 0.5],
// }
////////////////////////////////////////////////

////////////////////////////////////////////////
//                   PMTS
////////////////////////////////////////////////

////////////////////////////////////////////////
// // FAKE PMTS (acrylic cubes)
// {
//   name: "GEO",
//   index: "fakepmt",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "inner",
//   type: "box",
//   position: [-90.0, 0.0, 115.0],
//   size: [70.5,70.5,14.5],
//   material: "acrylic_berkeley", //quartz, acrylic_berkeley
//   color: [0.1, 0.3, 0.8, 0.1],
// }
/////////////////////////////////////////////////

///////////////////////////////
// COSMIC TAGS
///////////////////////////////

{
  name: "GEO",
  index: "tag1",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "tube",
  position: [-400.0, -400.0, -112.5],
  r_max: 5.0,
  size_z: 25.0,
  material: "air",
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
  position: [-400.0, -400.0, -257.5],
  r_max: 5.0,
  size_z: 25.0,
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}
/////////////////////////////////


/////////////////////////////////////////////////
// // ANALYSIS TUBES

//Container for pmts
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
  index: "small_pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
//  mother: "vessel",
  type: "pmtarray",
  pmt_model: "h11934", //h11934
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_CROSS_SIDE",
  orientation: "manual",
}

{
  name: "GEO",
  index: "big_pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
  type: "pmtarray",
  pmt_model: "r7081_hqe", //r7081_hqe, r11780_hqe
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_CLOSE",
  orientation: "point",
  orient_point: [-400.0, -400.0, -200.0],
}

// {
//   name: "GEO",
//   index: "trigger",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   mother: "inner",
//   type: "pmtarray",
//   pmt_model: "fast_test", //r7081_hqe, r11780_hqe, fast_test
//   pmt_detector_type: "idpmt",
//   sensitive_detector: "/mydet/pmt/inner",
//   efficiency_correction: 1.0,
//   pos_table: "PMTINFO_TRIGGER",
//   orientation: "point",
//   orient_point: [0.0, 0.0, 400.0],
// }

///////////////////////////////////////////////////////
////SOURCES
//Cosmics
//{
//  name: "GEO",
//  index: "cosmics",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "inner",
//  type: "tube",
//  r_max: 10.0,
//  size_z: 1.0,
//  position: [0.0, 0.0, 280.0],
//  material: "air",
//  color: [1.0, 0.1, 0.5, 0.5],
//}

//Radiactive
// {
//   name: "GEO",
//   index: "source",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "inner",
//   type: "tube",
//   position: [0.0, 0.0, 423.0],
//   r_max: 12.7,
//   size_z: 1.5, //half-height
//   material: "strontium", //strontium
//   color: [0.1, 1.0, 0.3, 0.8],
// }
