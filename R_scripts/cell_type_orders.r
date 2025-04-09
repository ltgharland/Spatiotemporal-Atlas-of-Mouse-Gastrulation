## Cell Type Lists Split by Germ Layer ##

ordered_refined_annotations <- c(
 #  "failed_QC",
    
  # Ectoderm
  "ExE ectoderm proximal",
  "ExE ectoderm distal",
#  "proximal ExE ectoderm #1",
#  "distal ExE ectoderm #2",
    
  #"ExE ectoderm #3",
  #"ExE ectoderm",
  "Otic placode",
  "Ectoderm",
  "Optic vesicle",
  "Hindbrain floor plate",
  "Hindbrain neural progenitors",
  "Ventral forebrain progenitors",
  #"Early dorsal forebrain progenitors",
  #"Late dorsal forebrain progenitors",
  "Midbrain/Hindbrain boundary",
  "Midbrain progenitors",
  "Dorsal midbrain neurons",
  "Ventral hindbrain progenitors",
  "Dorsal hindbrain progenitors",
    
  "Neural tube",
  "Migratory neural crest",
      "Spinal cord progenitors",
  "Dorsal spinal cord progenitors",
  "Branchial arch neural crest",
  "Frontonasal mesenchyme",

  "Non-neural ectoderm",
  "Surface ectoderm",
  "Epidermis",
  "Limb ectoderm",
  "Amniotic ectoderm",
  "Placodal ectoderm",
  "Otic neural progenitors",

  # Mesoderm
  "PGC",
  "Epiblast",
  "Rostral ectoderm",
  "Caudal epiblast",
  "Primitive Streak",
  "NMPs/Mesoderm-biased",
  "NMPs",
  "Nascent mesoderm",
   
# Cranial
     
  "Cranial mesoderm",
      "Embryo proper mesothelium",
    
# LPM 
  "Intermediate mesoderm",
 "ExE mesoderm",
 "ExE mesoderm and Anterior LPM",
    "Anterior LPM",
    "ExE mesoderm and Posterior LPM",
     # "allantois",
    #   "splanchnic LP mesoderm #1",
 # "splanchnic LP mesoderm #2",
   "Limb mesoderm",
 #  "Lateral plate mesoderm",
  #    "Forelimb",
  "Kidney primordium", 
    "Allantois",
     "Mesenchyme",
    
# Cardiac
          "Epicardium",
      "Pharyngeal mesoderm",
#  "Cardiopharyngeal progenitors",
  "Anterior cardiopharyngeal progenitors",
    
  "Cardiopharyngeal progenitors FHF",
#  "Cardiopharyngeal progenitors SHF",
  "Cardiomyocytes FHF 1",
  "Cardiomyocytes FHF 2",

  "Cardiomyocytes SHF 1",
  "Cardiomyocytes SHF 2",

        # Blood
      "Haematoendothelial progenitors",

  "Blood progenitors",
          "Chorioallantoic-derived erythroid progenitors",
  "Erythroid",

  #"Megakaryocyte progenitors",
  "MEP",
  "EMP",

    # Endothelium
  "Embryo proper endothelium #1",
  "Embryo proper endothelium #2",
  "Venous endothelium",
  "Endocardium #1",
  "Endocardium #2",
    "Allantois endothelium",


    
# Somitic
      "Caudal mesoderm",
     "Presomitic mesoderm",
  "Somitic mesoderm",
  "Posterior somitic tissues",
  "Paraxial mesoderm",
      "Anterior somitic tissues",
    "Endotome",
  "Sclerotome",
  "Dermomyotome",

    #YS
  "YS mesothelium",
  "YS endothelium",
  "YS mesothelium-derived endothelial progenitors",

      # Endoderm
  "Node",
  "Notochord",
  "Anterior Primitive Streak",
  "Gut tube",
  "Hindgut",
  #"Midgut",
  "Foregut",
  "Pharyngeal endoderm",
  "Thyroid primordium",
      "ExE endoderm",
  #"Parietal endoderm",
  "Visceral endoderm"
)

ordered_refined_celltypes_figure_2 <- c(
  'Epiblast', 'Rostral ectoderm', 'Primitive Streak', 'Anterior Primitive Streak', 
  'Nascent mesoderm', 'ExE ectoderm distal', 'ExE ectoderm proximal', 'Visceral endoderm', 
  'ExE endoderm', 'Node', 'PGC', 'Gut tube', 'ExE mesoderm', 
  'Haematoendothelial progenitors', 'Blood progenitors', 'ExE mesoderm and Anterior LPM', 
  'Non-neural ectoderm', 'Ectoderm', 'Paraxial mesoderm', 'Intermediate mesoderm', 
  'Caudal epiblast', 'Mesenchyme', 'Cardiopharyngeal progenitors FHF', 'Pharyngeal mesoderm', 
  'Anterior cardiopharyngeal progenitors', 'Cardiomyocytes FHF 1', 'Cardiomyocytes FHF 2', 
  'Cardiomyocytes SHF 1', 'Epicardium', 'ExE mesoderm and Posterior LPM', 'Anterior LPM', 
  'Limb mesoderm', 'Kidney primordium', 'Allantois', 'Allantois endothelium', 
  'Embryo proper endothelium #1', 'Embryo proper endothelium #2', 'Endocardium #1', 
  'Endocardium #2', 'Erythroid', 'Chorioallantoic-derived erythroid progenitors', 'EMP', 
  'MEP', 'Presomitic mesoderm', 'Somitic mesoderm', 'Anterior somitic tissues', 
  'Posterior somitic tissues', 'Dermomyotome', 'Endotome', 'Sclerotome', 
  'Embryo proper mesothelium', 'NMPs', 'Cranial mesoderm', 'Frontonasal mesenchyme', 
  'Foregut', 'Hindgut', 'Pharyngeal endoderm', 'Thyroid primordium', 'Notochord', 
  'Amniotic ectoderm', 'Surface ectoderm', 'Epidermis', 'Placodal ectoderm', 
  'Otic placode', 'Otic neural progenitors', 'Limb ectoderm', 'Migratory neural crest', 
  'Branchial arch neural crest', 'Neural tube',
  'Optic vesicle', 'Ventral forebrain progenitors', 'Midbrain progenitors', 
  'Dorsal midbrain neurons', 'Midbrain/Hindbrain boundary', 'Dorsal hindbrain progenitors', 
  'Hindbrain floor plate', 'Hindbrain neural progenitors', 'Ventral hindbrain progenitors', 
  'Dorsal spinal cord progenitors', 'Spinal cord progenitors'
)

Harland_seqFISH_early_order <- c(
  # Ectoderm
  "low quality",
  "ExE Ectoderm Proximal #1",
  "ExE Ectoderm Proximal #2",
  "ExE Ectoderm Distal #1",
  "ExE Ectoderm Distal #2",
    
  # EarlyMesoderm
  "Epiblast",
  "Rostral Ectoderm",
  "Caudal Epiblast",
  "Primitive Streak",
  #"anterior primitive streak",
  "Axial Mesendoderm",
  "Definitive Endoderm",
  "Mesodermal Wings",
  "ExE Mesoderm",
    
  # Blood
  "Haematoendothelial Progenitors #1",
  "Haematoendothelial Progenitors #2",

  # Endoderm
  "ExE Endoderm",
  "Visceral Endoderm #1",
  "Visceral Endoderm #2"
)

Harland_seqFISH_early_order_2 <- c(
  # Early Mesoderm
  "Epiblast",
  "Rostral Ectoderm",
  "Caudal Epiblast",
  "Primitive Streak",
  "Mesodermal Wings",
  "ExE Mesoderm",

  # Blood
  "Haematoendothelial Progenitors #1",
  "Haematoendothelial Progenitors #2",

  # Endoderm
  "Axial Mesendoderm",
  "Definitive Endoderm",
  "ExE Endoderm",
  "Visceral Endoderm #1",
  "Visceral Endoderm #2",

  # Ectoderm
  "ExE Ectoderm Proximal #1",
  "ExE Ectoderm Proximal #2",
  "ExE Ectoderm Distal #1",
  "ExE Ectoderm Distal #2",

  # Low Quality
  "low quality"
)