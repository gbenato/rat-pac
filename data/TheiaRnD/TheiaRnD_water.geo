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
  position: [-400.0, -400.0, -200.0],
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
  material: "chsrc_uvt_acrylic", //mirror
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
  material: "acrylic_black", //acrylic_black
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
  material: "acrylic_black", //acrylic_black
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
//  VESSEL
/////////////////
{
  name: "GEO",
  index: "container",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  invisible: 0, // omitted for visualization
  mother: "inner",
  type: "tube",
  r_max: 50.0,
  size_z: 15.0,
  position: [-400.0, -400.0, -152.5],
  material: "chsrc_uvt_acrylic",
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
  r_max: 45.0, //LAB->35.0, WBLS->30.0
  size_z: 14.0,
  position: [0.0, 0.0, 1.0],
  material: "water", //water, wbls_1pct_mod, wbls_1pct, wbls_5pct, wbls_10pct, labppo_scintillator
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
  pmt_model: "h11934", //h11934, r7600u-20
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
  pmt_model: "r7081_hqe", //r7081_hqe, r11780_hqe
  pmt_detector_type: "idpmt",
  sensitive_detector: "/mydet/pmt/inner",
  efficiency_correction: 1.0,
  pos_table: "PMTINFO_CLOSE",
  orientation: "point",
  orient_point: [-400.0, -400.0, -200.0],
}

{
  name: "GEO",
  index: "muon_pmts",
  valid_begin: [0, 0],
  valid_end: [0, 0],
  mother: "world",
  type: "pmtarray",
  pmt_model: "h11934", //h11934, r7081_hqe, r11780_hqe
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
   pmt_model: "h11934", //h11934, r7081_hqe, r11780_hqe, fast_test
   pmt_detector_type: "idpmt",
   sensitive_detector: "/mydet/pmt/inner",
   efficiency_correction: 1.0,
   pos_table: "PMTINFO_TRIGGER",
   orientation: "point",
   orient_point: [-400.0, -400.0, -152.5],
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
