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
//  size: [400.0,400.0,250.0], //mm, half-lenght
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
//  size: [370.0,370.0,220.0], //mm, half-lenght
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}

///////////////////////////////////////
//             VESSELS               //
///////////////////////////////////////

///////////////////////////////////////
//// FISH TANK
{
  name: "GEO",
  index: "outer_vessel",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  //Cylinder
  //  type: "tube",
  //  r_max: 150.0, //mm 200.0
  //  size_z: 50.0, //mm 19.0
  //Cuboid
  type: "box",
  size: [160.0,160.0,32.5],
  position: [0.0, 0.0, 0.0], //179.5
  material: "acrylic_berkeley", //quartz, acrylic_berkeley
  color: [0.1, 0.3, 0.8, 0.1],
}

// {
//   name: "GEO",
//   index: "inner_vessel",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "outer_vessel",
// //Cylinder
// //  type: "tube",
// //  r_max: 140.0, //mm
// //  size_z: 44.5, //mm
// //Cuboid
//   type: "box",
//   size: [190.0,190.0,20.0], //tank
//   position: [0.0, 0.0, 4.5],
//   material: "acrylic_berkeley", //water, acrylic_berkeley
//   color: [0.1, 0.8, 0.8, 0.1],
// }

{
  name: "GEO",
  index: "hollow_inner",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "outer_vessel",
  type: "tube",
  position: [0.0, 0.0, 0.0], //-15.0
  r_max: 5.0,
  size_z: 32.5,
  material: "air", //air
  color: [0.0, 0.0, 0.0, 0.1],
}

// {
//   name: "GEO",
//   index: "hollow_outer",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "outer_vessel",
//   type: "tube",
//   position: [0.0, 0.0, -45.0], //(-20.0)
//   r_max: 10.0, //5.0
//   size_z: 5.0, //5.0
//   material: "air", //air
//   color: [0.0, 0.0, 0.0, 0.1],
// }

{
  name: "GEO",
  index: "container",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner", //inner_vessel
  //Cylinder
  //  type: "tube",
  //  r_max: 30.0, //mm
  //  size_z: 15.0, //mm 19.0
  //Cuboid
  type: "box",
  rotation:  [0.0, 0.0, 0.0],
  size: [100.0,100.0,15.0], //17.5
  position: [0.0, 0.0, 47.5],
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
  //Cylinder
  type: "tube",
  r_max: 50.0, //10.0, 30.0
  size_z: 13.0, //mm 17.5
  //Cuboid
  //  type: "box",
  //  size: [25.0,25.0,10.0], //25.0, 10.0
  position: [0.0, 0.0, 0.0], //2.5
  material: "labppo_scintillator", //water, wbls_5pct, wbls_10pct, labppo_scintillator
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

/////////////////////////////////////////////////
// ANALYSIS TUBES

//Container for pmts
//{
//  name: "GEO",
//  index: "pmt_holder",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "inner",
////Cuboid
//  type: "box",
//  size: [121.0,121.0,15.0],
//  position: [0.0, 0.0, -40.0],
//  material: "acrylic_black",
//  color: [0.5, 0.0, 0.0, 0.5],
//}

{
  name: "GEO",
  index: "pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
//  mother: "vessel",
  type: "pmtarray",
  pmt_model: "h11934", //h11934, r7081_hqe, r11780_hqe
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_ARRAY",
  orientation: "manual",
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

///////////////////////////////////////////////////////
////SOURCES
//Cosmics
{
  name: "GEO",
  index: "cosmics",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "tube",
  r_max: 10.0,
  size_z: 1.0,
  position: [0.0, 0.0, 280.0],
  material: "air",
  color: [1.0, 0.1, 0.5, 0.5],
}
