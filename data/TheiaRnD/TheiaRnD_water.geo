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
  material: "chsrc_uvt_acrylic",
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
  position: [0.0, 0.0, 0.5],
  r_max: 5.0,
  size_z: 31.0,
  material: "acrylic_black",
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
  size_z: 1.0,
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}

/////////////////
//  VESSEL
/////////////////
{
  name: "GEO",
  index: "vessel",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "CheSSVessel",
  r_max: 50.0,
  size_z: 15.0,
  rotation:  [0.0, 0.0, 135.0],
  position: [-398.0, -367.0, -220.91],
  material: "chsrc_uvt_acrylic",
  color: [0.1, 0.3, 0.8, 0.1],
}

{
  name: "GEO",
  index: "cavity",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "vessel",
  type: "tube",
  r_max: 20.0,
  size_z: 1.905,
  position: [0.0, 0.0, 14.685],
  material: "air",
  color: [0.0, 0.0, 0.0, 0.1],
}

{
  name: "GEO",
  index: "content",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "vessel",
  type: "tube",
  r_max: 45.0, //LAB->45mm, WBLS->35mm
  size_z: 12.89,
  position: [0.0, 0.0, -1.7],
  material: "water",
  color: [0.5, 0.1, 0.5, 0.5],
}

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
  position: [-398.0, -367.0, -182.987],
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
  pmt_model: "h11934",
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
  pmt_model: "r7081_hqe",
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_CLOSE",
  orientation: "point",
  orient_point: [-398.0, -367.0, -237.0],
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
   index: "trigger_pmt",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   mother: "inner",
   type: "pmtarray",
   pmt_model: "h11934",
   pmt_detector_type: "idpmt",
   sensitive_detector: "/mydet/pmt/inner",
   efficiency_correction: 1.0,
   pos_table: "PMTINFO_TRIGGER",
   orientation: "point",
   orient_point: [-398.0, -367.0, -222.2],
}

{
   name: "GEO",
   index: "panels_pmt",
   valid_begin: [0, 0],
   valid_end: [0, 0],
   mother: "world",
   type: "pmtarray",
   pmt_model: "h11934",
   pmt_detector_type: "idpmt",
   sensitive_detector: "/mydet/pmt/inner",
   efficiency_correction: 1.0,
   pos_table: "PMTINFO_PANELS",
   orientation: "manual",
   orient_point: [0.0, 0.0, 400.0],
}
