{
  name: "GEO",
  index: "world",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 1, // omitted for visualization
  mother: "",
  type: "box",
  size: [3000.0,3000.0,3000.0], //mm, half-lenght
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

{
  name: "GEO",
  index: "pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
  type: "pmtarray",
  pmt_model: "r7081_hqe", //r7081_hqe, r11780_hqe, fast_test
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_CLOSE",
  orientation: "point",
  orient_point: [-400.0, -400.0, 0.0],
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
//   efficiency_correction: 1.027,
//   pos_table: "PMTINFO_TRIGGER",
//   orientation: "point",
//   orient_point: [0.0, 0.0, 400.0],
// }

{
  name: "GEO",
  index: "vessel",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  position: [-400.0, -400.0, 0.0],
  size: [12.7,25.4,50.8], //mm, half-lenght
  rotation:  [0.0, 0.0, -45.0],
  material: "acrylic_berkeley", //acrylic_berkeley
  color: [0.1, 0.3, 0.8, 0.1],
}

//{
//  name: "GEO",
//  index: "content",
//  valid_begin: [0, 0],
//  valid_end: [0, 0],
//  invisible: 0, // omitted for visualization
//  mother: "vessel",
//  type: "tube",
//  position: [0.0, 0.0, 10.0],
//  r_max: 10.0,
//  size_z: 40.8,
//  material: "scintillator", //air, water, scintillator
//  color: [0.1, 0.3, 0.8, 0.3],
//}

{
  name: "GEO",
  index: "envelope",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "tube",
  position: [-410.0, -410.0, 0.0],
  rotation:  [-90.0, 45.0, 0.0],
  r_max: 12.808,
  size_z: 1.587, //half height
  material: "acrylic_black", //strontium
  color: [0.1, 1.0, 0.3, 0.8],
}

{
  name: "GEO",
  index: "source",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "envelope",
  type: "tube",
  position: [0.0, 0.0, 0.254],
  r_max: 3.175,
  size_z: 1.333, //half height
  material: "strontium", //strontium
  color: [0.1, 1.0, 1.0, 0.8],
}

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


// {
//   name: "GEO",
//   index: "ruler",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "inner",
//   type: "sphere",
//   position: [-400.0, -400.0, 0.0],
//   r_max: 290.0,
//   material: "air",
//   color: [0.1, 0.3, 0.8, 0.1],
// }
