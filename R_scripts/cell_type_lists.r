## Cell Type Lists Split by Germ Layer ##

mesoderm_layer <- c(
  "Pharyngeal mesoderm",
  "Cardiopharyngeal progenitors",
  "Cardiopharyngeal progenitors SHF",
  "Anterior cardiopharyngeal progenitors",
  "Cardiopharyngeal progenitors FHF",
  "Cardiomyocytes FHF 1",
  "Cardiomyocytes FHF 2",
  "Cardiomyocytes SHF 1",
  "Cardiomyocytes SHF 2",
  "Epicardium",
  "Endocardium",
  "Intermediate mesoderm",
  "ExE mesoderm and Anterior LPM",
  "ExE mesoderm and Posterior LPM",
  "Anterior LPM",
  "Limb mesoderm",
  "Kidney primordium",
  "Allantois",
  "Paraxial mesoderm",
  "Presomitic mesoderm",
  "Somitic mesoderm",
  "Anterior somitic tissues",
  "Posterior somitic tissues",
  "Endotome",
  "Dermomyotome",
  "Cranial mesoderm and anterior somitic tissue",
  "Sclerotome",
  "Embryo proper mesothelium",
  "Mesenchyme",
  "NMPs",
  "Cranial mesoderm",
  "Frontonasal mesenchyme",
  "Allantois endothelium",
  "Embryo proper endothelium #1",
  "Embryo proper endothelium #2",
  "Endocardium #1",
  "Endocardium #2",
  "Haematoendothelial progenitors",
  "Blood progenitors",
  "Erythroid",
  "Chorioallantoic-derived erythroid progenitors",
  "EMP",
  "MEP"
)

ectoderm_layer <- c(
  "Amniotic ectoderm",
  "Surface ectoderm",
  "Epidermis",
  "Placodal ectoderm",
  "Otic placode",
  "Otic neural progenitors",
  "Limb ectoderm",
  "Migratory neural crest",
  "Branchial arch neural crest",
  "Ectoderm",
  "Non-neural ectoderm",
  "Neural tube",
  "Optic vesicle",
  "Early dorsal forebrain progenitors",
  "Late dorsal forebrain progenitors",
  "Ventral forebrain progenitors",
  "Midbrain progenitors",
  "Dorsal midbrain neurons",
  "Midbrain/Hindbrain boundary",
  "Dorsal hindbrain progenitors",
  "Hindbrain floor plate",
  "Hindbrain neural progenitors",
  "Ventral hindbrain progenitors",
  "Dorsal spinal cord progenitors",
  "Spinal cord progenitors",
  "Ventral neural tube"
)

somitic_layer <- c(
  "Anterior somitic tissues",
  "Posterior somitic tissues",
  "Dermomyotome",
  "Sclerotome",
  "Endotome",
  "Somitic mesoderm"
)

endoderm_layer <- c(
  # Primitive Endoderm
  "Parietal endoderm",
  "ExE endoderm",
  "Visceral endoderm",
  
  # Definitive Endoderm
  "Anterior Primitive Streak",
  "Gut tube",
  "Foregut",
  "Midgut",
  "Hindgut",
  "Pharyngeal endoderm",
  "Thyroid primordium",
  
  # Additional Structures
  "Node",
  "Notochord"
)

## Cell Type Lists Split by Tissue ##

## E8.5 Cell Types ##

LPM <- c(
  "Pharyngeal mesoderm",
  "Anterior cardiopharyngeal progenitors",
  "Cardiopharyngeal progenitors FHF",
  "Cardiomyocytes FHF 1",
  "Cardiomyocytes FHF 2",
  "Cardiomyocytes SHF 1",
  "Epicardium",
  "Intermediate mesoderm",
  "ExE mesoderm and Anterior LPM",
  "ExE mesoderm and Posterior LPM",
  "Anterior LPM",
  "Limb mesoderm",
  "Kidney primordium",
  "Allantois"
)

cranial_somitic <- c(
  "Paraxial mesoderm",
  "Presomitic mesoderm",
  "Somitic mesoderm",
  "Anterior somitic tissues",
  "Posterior somitic tissues",
  "Dermomyotome",
  "Endotome",
  "Sclerotome",
  "Embryo proper mesothelium",
  "Mesenchyme",
  "NMPs",
  "Cranial mesoderm",
  "Frontonasal mesenchyme"
)

hematovascular <- c(
  "Allantois endothelium",
  "Embryo proper endothelium #1",
  "Embryo proper endothelium #2",
  "Endocardium #1",
  "Endocardium #2",
  "Haematoendothelial progenitors",
  "Blood progenitors",
  "Erythroid",
  "Chorioallantoic-derived erythroid progenitors",
  "EMP",
  "MEP"
)

endoderm <- c(
  "Anterior Primitive Streak",
  "Gut tube",
  "Foregut",
  "Hindgut",
  "Pharyngeal endoderm",
  "Thyroid primordium",
  "ExE endoderm",
  "Visceral endoderm",
  "Notochord",
  "Node"
)

surface_ectoderm <- c(
  "Amniotic ectoderm",
  "Surface ectoderm",
  "Epidermis",
  "Placodal ectoderm",
  "Otic placode",
  "Otic neural progenitors",
  "Limb ectoderm",
  "Migratory neural crest",
  "Branchial arch neural crest",
  "Ectoderm",
  "Non-neural ectoderm"
)

neural_tube <- c(
  "Neural tube",
  "Optic vesicle",
  "Ventral forebrain progenitors",
  "Midbrain progenitors",
  "Dorsal midbrain neurons",
  "Midbrain/Hindbrain boundary",
  "Dorsal hindbrain progenitors",
  "Hindbrain floor plate",
  "Hindbrain neural progenitors",
  "Ventral hindbrain progenitors",
  "Dorsal spinal cord progenitors",
  "Spinal cord progenitors",
  "Ventral neural tube"
)

## E7.25 and E7.5 Cell Types ##

early_gastrula <- c(
  "Anterior Primitive Streak",
  "Parietal endoderm",
  "Gut tube",
  "ExE endoderm",
  "Visceral endoderm",
  "ExE ectoderm distal",
  "ExE ectoderm proximal",
  "Epiblast",
  "Primitive Streak",
  "Nascent mesoderm",
  "PGC",
  "ExE mesoderm",
  "Haematoendothelial progenitors",
  "Blood progenitors",
  "ExE mesoderm and Anterior LPM",
  "Non-neural ectoderm",
  "Ectoderm",
  "Intermediate mesoderm",
  "Caudal epiblast",
  "Rostral ectoderm"
)