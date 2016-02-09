///////////////////////////////
//
// PMT configuration with only one big PMT close from a piece of acrylic
// and one trigger PMT attached to the acrylic. A source of strontium is
// attached at its back.
//
///////////////////////////////


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
  index: "base",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "box",
  position: [-450.0, 0.0, -388.9],
  size: [63.5,114.3,63.5], //mm, half-lenght
  material: "styrofoam", //styrofoam
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
  position: [-450.0, 0.0, -300.0],
  size: [12.7,50.8,25.4], //mm
//  size: [2.5,50.8,25.4], //mm
  material: "acrylic_berkeley", //labppo_scintillator, scintillator, acrylic_berkeley
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
//   position: [-392.8, 0.0, -300.0],
//   size: [44.5,10.,10.0], //mm, half-lenght
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
  position: [-464.287, -25.0, -300.0], //-483.587
  rotation:  [0.0, 90.0, 0.0],
  r_max: 12.808,
  size_z: 1.587, //half height
  material: "acrylic_black",
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

{
  name: "GEO",
  index: "pmt_trig",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
  type: "pmtarray",
  pmt_model: "h11934",//h11934
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_TRIGGER",
  orientation: "manual",
  orient_point: [-464.287, -25.0, -300.0],
}

{
  name: "GEO",
  index: "pmt",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "inner",
  type: "pmtarray",
  pmt_model: "r7081_hqe", //r7081_hqe, h11934
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_ONEPMT",
  orientation: "manual",
  orient_point: [-450.0, 0.0, -300.0],
}
