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
  size: [500.0,250.0,250.0], //mm, half-lenght
  material: "acrylic_black", //cardboard
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
  size: [480.0,230.0,230.0], //mm, half-lenght
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "base",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  position: [-225.0, 0.0, -75.0],
  size: [50.0,150.0,50.0], //mm, half-lenght
  material: "styrofoam",
  color: [1.0, 1.0, 1.0, 1.0],
}

{
  name: "GEO",
  index: "vessel",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  position: [-180.0, 0.0, 0.0],
  rotation:  [45.0, 0.0, 0.0],
  size: [2.5,25.4,25.4], //mm, half-lenght
  material: "acrylic_berkeley",
  color: [0.1, 0.3, 0.8, 0.1],
}

// {
//   name: "GEO",
//   index: "ruler",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   invisible: 0, // omitted for visualization
//   mother: "inner",
//   type: "box",
//   position: [-172.75, 0.0, 0.0],
//   size: [4.75,30.,30.0], //mm, half-lenght
//   material: "air", //acrylic_uvt
//   color: [0.1, 0.3, 0.8, 0.1],
// }

{
  name: "GEO",
  index: "envelope",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "tube",
  position: [-183.5, 0.0, 0.0],
  rotation:  [0.0, 90.0, 0.0],
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
//   index: "pmt_trig",
//   valid_begin: [0, 0],
//   valid_end: [0, 0],
//   mother: "inner",
//   type: "pmtarray",
//   pmt_model: "fast_test",// r7081_hqe, r11780_hqe
//   pmt_detector_type: "idpmt",
//   sensitive_detector: "/mydet/pmt/inner",
//   efficiency_correction: 1.0,
//   pos_table: "PMTINFO_TRIGGER",
//   orientation: "point",
//   orient_point: [-180.0, 0.0, 0.0],
// }

{
  name: "GEO",
  index: "pmt",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
  type: "pmtarray",
  pmt_model: "r7081_hqe",// r7081_hqe, r11780_hqe
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 0.600000,
  pos_table: "PMTINFO",
  orientation: "point",
  orient_point: [-180.0, 0.0, 0.0],
}
