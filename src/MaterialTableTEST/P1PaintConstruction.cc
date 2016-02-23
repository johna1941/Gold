// Test to see if I can write code for the paint properties.


// These need to be subbed out for the true values.
  G4double photonEnergy2[] =
  { 1*eV, 1*eV, 1*eV, 1*eV,
    1*eV, 1*eV, 1*eV, 1*eV,
    1*eV, 1*eV, 1*eV, 1*eV,
    1*eV, 1*eV, 1*eV, 1*eV,
    1*eV, 1*eV, 1*eV, 1*eV,
    1*eV, 1*eV, 1*eV, 1*eV,
    1*eV, 1*eV, 1*eV, 1*eV,
    1*eV, 1*eV, 1*eV, 1*eV };

  G4double RIndexPaint[] =
  { 1, 1, 1, 1, 
    1, 1, 1, 1, 
    1, 1, 1, 1, 
    1, 1, 1, 1, 
    1, 1, 1, 1,
    1, 1, 1, 1, 
    1, 1, 1, 1, 
    1, 1, 1, 1};

  G4double AbsPaint[] =
 { 1*m, 1*m, 1*m, 1*m,
   1*m, 1*m, 1*m, 1*m, 
   1*m, 1*m, 1*m, 1*m,
   1*m, 1*m, 1*m, 1*m, 
   1*m, 1*m, 1*m, 1*m, 
   1*m, 1*m, 1*m, 1*m,
   1*m, 1*m, 1*m, 1*m, 
   1*m, 1*m, 1*m, 1*m };

   // Paint is not a scintillator, this should go:
   /*
  G4double scintilFastPaint[] =
  { 1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00 };
  G4double scintilSlowPaint[] =
  { 1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 
    1.00, 1.00, 1.00, 1.00 };
    */
  // Health check
  const G4int nEntries2 = sizeof(photonEnergy2)/sizeof(G4double);
  assert(sizeof(RIndexPaint) == sizeof(photonEnergy));
  assert(sizeof(AbsPaint) == sizeof(photonEnergy));
  // assert(sizeof(scintilFastPaint) == sizeof(photonEnergy));
  // assert(sizeof(scintilSlowPaint) == sizeof(photonEnergy));
// Create material properties table and add properties
  G4MaterialPropertiesTable* paint_mpt = new G4MaterialPropertiesTable();
  // Add to material properties table
  paint_mpt->AddProperty("RINDEXPAINT",       photonEnergy2, RIndexPaint,   nEntries)
  ->SetSpline(true);
  paint_mpt->AddProperty("ABSLENGTHPAINT",    photonEnergy2, AbsPaint,      nEntries)
  ->SetSpline(true);
  /*
  paint_mpt->AddProperty("FASTCOMPONENTPAINT",photonEnergy2, scintilFastPaint,     nEntries)
  ->SetSpline(true);
  paint_mpt->AddProperty("SLOWCOMPONENTPAINT",photonEnergy2, scintilSlowPaint,     nEntries)
  ->SetSpline(true);
  */

  // Scintillator constant properties:
  /*
  paint_mpt->AddConstProperty("SCINTILLATIONYIELD",0.1/MeV);
  paint_mpt->AddConstProperty("RESOLUTIONSCALE",1.0);
  paint_mpt->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  paint_mpt->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  paint_mpt->AddConstProperty("YIELDRATIO",0.8);
  */

  G4cout << "Paint G4MaterialPropertiesTable\n"; paint_mpt->DumpTable();
  // Associate material properties table with the liquid scintillator material
  paint->SetMaterialPropertiesTable(paint_mpt);

  // Optical properties of the surface of the scintillator
  G4OpticalSurface* paint_surface = new G4OpticalSurface("paint_surface");
  paint_surface->SetType(dielectric_dielectric); // If both surfaces have refractive properties added, this will actually calculate reflection for us
  paint_surface->SetFinish(groundfrontpainted);
  paint_surface->SetModel(unified);
  G4cout << "paint_surface\n"; paint_surface->DumpInfo();
  // Create material properties table and add properties
  if (fReflectivity < 0.) {
    G4cout << "Reflectivity not set!" << G4endl;
    abort();
  }

  // Skin of variable reflective for the messenger, I think
  G4double reflectivity[nEntries]; for (auto& r: reflectivity) r = fReflectivity;
  G4MaterialPropertiesTable* mptForSkin = new G4MaterialPropertiesTable();  
  mptForSkin->AddProperty("REFLECTIVITY", photonEnergy, reflectivity, nEntries)
  ->SetSpline(true);
  G4cout << "Skin G4MaterialPropertiesTable\n"; mptForSkin->DumpTable();
  // Associates the material properties with the surface of the liquid scintillator. 
  scint_surface->SetMaterialPropertiesTable(mptForSkin); 
