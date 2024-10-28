#######################################
# Muon Collider Detector TARGET model
#
# Michele Selvaggi michele.selvaggi@cern.ch
# Ulrike Schnoor ulrike.schnoor@cern.ch
#
#
# !!! DISCLAIMER !!!
#
# The parameterisation of the Muon Collider
# has to be intended as a target performance.
# This has not been validated by full simulation.
# Hybrid between FCC-hh and CLIC performance.
#
#
#######################################

# Magnetic field in the solenoid
set B 5.0
#radius of solenoid CHECK!!!!
set R 2.055
#half-length of solenoid
set HL 2.509

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
    ParticlePropagator

    ChargedHadronTrackingEfficiency
    ElectronTrackingEfficiency
    MuonTrackingEfficiency

    TrackMergerPre

    TrackSmearing

    TrackMerger

    ECal
    HCal

    Calorimeter
    EFlowMerger

    PhotonEfficiency
    PhotonIsolation

    ElectronFilter
    ElectronEfficiency
    ElectronIsolation

    ChargedHadronFilter
    MuonFilter

    MuonEfficiency
    MuonIsolation

    EFlowFilter

    NeutrinoFilter


    GenJetFinder
    FastJetFinderKt

    JetFlavorAssociation
    BTaggingWP70

    TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
    set InputArray Delphes/stableParticles

    set OutputArray stableParticles
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons
    set MuonOutputArray muons

    # radius of the magnetic field coverage in the calorimeter, in m
    set Radius $R
    # half-length of the magnetic field coverage in the calorimeter, in m
    set HalfLength $HL

    # magnetic field, in T
    set Bz $B
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray ParticlePropagator/chargedHadrons
    set OutputArray chargedHadrons
    # Current full simulation with CLICdet provides for pions:

    set EfficiencyFormula {
 (pt <= 0.1)                                                                        * (0.000) +
 (pt > 0.1)                                * (abs(eta) > 2.5)                      * (0.000) +
 (pt > 0.1) * (energy >= 80)               * (abs(eta) < 2.5)                      * (1.000) +
 (pt > 0.1) * (energy < 80 && energy >= 3) * (abs(eta) <=2.5 && abs(eta) > 2.34)   * (0.994) +
 (pt > 0.1) * (energy < 80 && energy >= 3) * (abs(eta) <= 2.34)                     * (1.000) +
 (pt > 0.1) * (energy < 3)                 * (abs(eta) <= 2.5 && abs(eta) > 0.55 ) * (0.990) +
 (pt > 0.1) * (energy < 3)                 * (abs(eta) <= 0.55 )                    * (1.000) 
    }
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
    set InputArray ParticlePropagator/electrons
    set OutputArray electrons


    # Current full simulation with CLICdet provides for electrons:
    set EfficiencyFormula {
 (pt <= 0.1)                                                                        * (0.000) +
 (pt > 0.1)                                * (abs(eta) > 2.5)                      * (0.000) +
 (pt > 0.1) * (energy >= 80)               * (abs(eta) <= 2.5 && abs(eta) > 2.44 ) * (0.993) +
 (pt > 0.1) * (energy >= 80)               * (abs(eta) <= 2.44 && abs(eta) > 2.34 ) * (0.997) +
 (pt > 0.1) * (energy >= 80)               * (abs(eta) <= 2.34  )                   * (1.000) +
 (pt > 0.1) * (energy < 80 && energy >= 5) * (abs(eta) <= 2.5 && abs(eta) > 2.17 ) * (0.998) +
 (pt > 0.1) * (energy < 80 && energy >= 5) * (abs(eta) <= 2.17)                     * (1.000) +
 (pt > 0.1) * (energy < 5)                 * (abs(eta) <= 2.5 && abs(eta) > 2.34 ) * (1.000) +
 (pt > 0.1) * (energy < 5)                 * (abs(eta) <= 2.34 && abs(eta) > 0.76 ) * (0.997) +
 (pt > 0.1) * (energy < 5)                 * (abs(eta) <= 0.76)                     * (0.999)
    }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
    set InputArray ParticlePropagator/muons
    set OutputArray muons

    # Current full simulation with CLICdet provides for muons:
    set EfficiencyFormula {
  (pt <= 0.1)                                                                          * (0.000) +
  (pt > 0.1) * (abs(eta) > 2.5)                                                       * (0.000) +
  (pt > 0.1) * (abs(eta) <= 2.5 && abs(eta) > 2.44 ) * (energy >= 80)                 * (0.994) +
  (pt > 0.1) * (abs(eta) <= 2.5 && abs(eta) > 2.44 ) * (energy >= 5 && energy < 80)   * (0.996) +
  (pt > 0.1) * (abs(eta) <= 2.5 && abs(eta) > 2.44 ) * (energy < 5 )                  * (0.996) +
  (pt > 0.1) * (abs(eta) <= 2.44 )                    * (energy >= 5 )                 * (1.000) +
  (pt > 0.1) * (abs(eta) <= 2.44 && abs(eta) > 2.25 ) * (energy < 5 )                  * (0.999) +
  (pt > 0.1) * (abs(eta) <= 2.25 )                    * (energy < 5 )                  * (1.000) 
    }
}

##############
# Track merger
##############

module Merger TrackMergerPre {
# add InputArray InputArray
  add InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  add InputArray ElectronTrackingEfficiency/electrons
  add InputArray MuonTrackingEfficiency/muons
  set OutputArray tracks
}



########################################
# Smearing for charged tracks
########################################

module TrackCovariance TrackSmearing {

    set InputArray TrackMergerPre/tracks
    set OutputArray tracks

    ## minimum number of hits to accept a track
    set NMinHits 6

    ## magnetic field
    set Bz $B

    ## scale factors
    #set ElectronScaleFactor  {1.0}


    set DetectorGeometry {
        1 PIPE -100 100 0.0228 0.0012 0.35276 0 0 0 0 0 0
        1 VTX -0.13 0.13 0.029 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        1 VTX -0.13 0.13 0.04 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        1 VTX -0.13 0.13 0.05 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        1 VTX -0.13 0.13 0.073 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        1 VTX -0.13 0.13 0.101 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        1 ITK -0.4816 0.4816 0.164 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        1 ITK -0.4816 0.4816 0.354 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        1 ITK -0.6923 0.6923 0.554 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        1 OTK -1.2642 1.2642 0.819 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        1 OTK -1.2642 1.2642 1.153 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        1 OTK -1.2642 1.2642 1.486 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        1 MAG -2.5 2.5 2.25 0.05 0.0658 0 0 0 0 0 0
        2 VTXDSK 0.065 0.112 -0.366 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 VTXDSK 0.053 0.112 -0.298 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 VTXDSK 0.041 0.112 -0.23 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 VTXDSK 0.032 0.112 -0.18 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 VTXDSK 0.032 0.112 0.18 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 VTXDSK 0.041 0.112 0.23 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 VTXDSK 0.053 0.112 0.298 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 VTXDSK 0.065 0.112 0.366 5e-05 0.0937 2 0 1.5708 5e-06 5e-06 1
        2 ITKDSK 0.277 0.56 -2.19 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.257 0.56 -1.946 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.239 0.56 -1.741 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.215 0.56 -1.457 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.191 0.56 -1.173 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.156 0.56 -0.888 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.106 0.41 -0.604 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.106 0.41 0.604 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.156 0.56 0.888 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.191 0.56 1.173 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.215 0.56 1.457 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.239 0.56 1.741 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.257 0.56 1.946 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 ITKDSK 0.277 0.56 2.19 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 -2.19 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 -1.933 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 -1.667 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 -1.41 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 1.41 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 1.667 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 1.933 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 2.19 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1
        2 OTKDSK 0.6175 1.4302 2.19 0.000956 0.0937 2 0 1.5708 7e-06 9e-05 1

    }

}

##############
# Track merger
##############

module Merger TrackMerger {

  add InputArray TrackSmearing/tracks
  set OutputArray tracks
}


#############
#   ECAL
#############

module SimpleCalorimeter ECal {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackMerger/tracks

    set TowerOutputArray ecalTowers
    set EFlowTrackOutputArray eflowTracks
    set EFlowTowerOutputArray eflowPhotons

    set IsEcal true

    set EnergyMin 0.5
    set EnergySignificanceMin 1.0

    set SmearTowerCenter true

    set pi [expr {acos(-1)}]

    # lists of the edges of each tower in eta and phi
    # each list starts with the lower edge of the first tower
    # the list ends with the higher edged of the last tower

    #ECAL barrel: dphi = 0.34 degree, deta=0.005 towers up to |eta| <=1.2
    #ECAL endcaps: dphi = 0.8 degree, deta=0.02 towers up to |eta| <=2.5

    #barrel:
    #dphi = 0.34 degree towers up to eta <=1.1 (1050 cells)
    set PhiBins {}
    for {set i -525} {$i <= 525} {incr i} {
	add PhiBins [expr {$i * $pi/525.0 }]
    }
    # 0.005 unit (10x10 mm^2) in eta up to eta <=1.1
    for {set i -221} {$i <=221} {incr i} {
	set eta [expr {$i * 0.005}]
	add EtaPhiBins $eta $PhiBins
    }

    #endcaps:
    #dphi = 0.92 degree towers for 1.2 < eta <=2.5
    set PhiBins {}
    for {set i -196} {$i <= 196} {incr i} {
	add PhiBins [expr {$i * $pi/196.}]
    }
    #deta=0.02 units for 1.2 < |eta| <=2.5
    #first, from -2.5 to -1.2, there will be (1.3/0.02=)65 segments
    for {set i 1} {$i <=166} {incr i} {
	set eta [expr {-2.42 + $i * 0.008}]
	add EtaPhiBins $eta $PhiBins
    }
    #same for 1.2 to 2.5
    for  {set i 1} {$i <=166} {incr i} {
	set eta [expr {1.092 + $i*0.008}]
	add EtaPhiBins $eta $PhiBins
    }


    # default energy fractions {abs(PDG code)} {fraction of energy deposited in ECAL}

    add EnergyFraction {0} {0.0}
    # energy fractions for e, gamma and pi0
    add EnergyFraction {11} {1.0}
    add EnergyFraction {22} {1.0}
    add EnergyFraction {111} {1.0}
    # energy fractions for muon, neutrinos and neutralinos
    add EnergyFraction {12} {0.0}
    add EnergyFraction {13} {0.0}
    add EnergyFraction {14} {0.0}
    add EnergyFraction {16} {0.0}
    add EnergyFraction {1000022} {0.0}
    add EnergyFraction {1000023} {0.0}
    add EnergyFraction {1000025} {0.0}
    add EnergyFraction {1000035} {0.0}
    add EnergyFraction {1000045} {0.0}
    # energy fractions for K0short and Lambda
    add EnergyFraction {310} {0.3}
    add EnergyFraction {3122} {0.3}

    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    # sqrt(energy^2*c^2 + energy*a^2) where a is the stochastic term and c the constant term
    set ResolutionFormula {
	(abs(eta) <= 2.5) * sqrt( energy^2*0.01^2 + energy*0.10^2  )}
}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray ECal/eflowTracks

    set TowerOutputArray hcalTowers
    set EFlowTrackOutputArray eflowTracks
    set EFlowTowerOutputArray eflowNeutralHadrons

    set IsEcal false

    set EnergyMin 1.0
    set EnergySignificanceMin 1.0

    set SmearTowerCenter true

    set pi [expr {acos(-1)}]

    # lists of the edges of each tower in eta and phi
    # each list starts with the lower edge of the first tower
    # the list ends with the higher edged of the last tower


    #HCAL barrel: dphi = 1 degree, deta= 0.02 towers up to |eta| <=0.8
    #HCAL ring: dphi = 1 degree, deta= 0.02 towers up to |eta| <=0.9
    #HCAL endcaps: dphi = 6 degree, deta = 0.1  up to |eta| <=2.5
    #HCAL cell sizes always 30x30 mm^2

    #barrel and ring:
    #dphi = 1 degree up to |eta| <=0.9
    set PhiBins {}
    for {set i -180} {$i <=180} {incr i} {
	add PhiBins [expr {$i * $pi/180.0}]
    }
    #deta= 0.02 towers up to |eta| <=0.9
    for {set i -45} {$i <=45} {incr i} {
	set eta [expr {$i * 0.02}]
	add EtaPhiBins $eta $PhiBins
    }

    #endcaps:
    # dphi = 6 degree
    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    # deta =0.1 for 0.9 < |eta| <=2.5
    #for -2.5 to -0.9, 21 segments
    for {set i 1} {$i <=17} {incr i} {
	set eta [expr {-2.5 + $i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }
    #same for 0.9 to 2.5
    for {set i 1} {$i <=17} {incr i} {
	set eta [expr {0.8 + $i * 0.1 }]
	add EtaPhiBins $eta $PhiBins
    }

    # default energy fractions {abs(PDG code)} {Fecal Fhcal}
    add EnergyFraction {0} {1.0}
    # energy fractions for e, gamma and pi0
    add EnergyFraction {11} {0.0}
    add EnergyFraction {22} {0.0}
    add EnergyFraction {111} {0.0}
    # energy fractions for muon, neutrinos and neutralinos
    add EnergyFraction {12} {0.0}
    add EnergyFraction {13} {0.0}
    add EnergyFraction {14} {0.0}
    add EnergyFraction {16} {0.0}
    add EnergyFraction {1000022} {0.0}
    add EnergyFraction {1000023} {0.0}
    add EnergyFraction {1000025} {0.0}
    add EnergyFraction {1000035} {0.0}
    add EnergyFraction {1000045} {0.0}
    # energy fractions for K0short and Lambda
    add EnergyFraction {310} {0.7}
    add EnergyFraction {3122} {0.7}

    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    # sqrt(energy^2*c^2 + energy*a^2) where a is the stochastic term and c the constant term
    set ResolutionFormula {
	    (abs(eta)<= 2.5) * sqrt(energy*0.308^2  + energy^2*0.050^2)
	}

}

#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
    set InputArray HCal/eflowTracks
    set OutputArray electrons
    set Invert true
    add PdgCode {11}
    add PdgCode {-11}
}

######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
    set InputArray HCal/eflowTracks
    set OutputArray chargedHadrons

    add PdgCode {11}
    add PdgCode {-11}
    add PdgCode {13}
    add PdgCode {-13}
}

#################
# Muon filter
#################

module PdgCodeFilter MuonFilter {
  set InputArray HCal/eflowTracks
  set OutputArray muons
  set Invert true
  add PdgCode {13}
  add PdgCode {-13}
}




###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger Calorimeter {
    # add InputArray InputArray
    add InputArray ECal/ecalTowers
    add InputArray HCal/hcalTowers
    set OutputArray towers
}


####################
# Energy flow merger
####################

module Merger EFlowMerger {
    # add InputArray InputArray
    add InputArray HCal/eflowTracks
    add InputArray ECal/eflowPhotons
    add InputArray HCal/eflowNeutralHadrons
    set OutputArray eflow
}

######################
# EFlowFilter
######################

# module PdgCodeFilter EFlowFilter {
#   set InputArray EFlowMerger/eflow
#   set OutputArray eflow

#   add PdgCode {11}
#   add PdgCode {-11}
#   add PdgCode {13}
#   add PdgCode {-13}
# }


###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
    set InputArray ECal/eflowPhotons
    set OutputArray photons

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # efficiency formula for photons
    set EfficiencyFormula {
	(energy < 2.0 ) * (0.000) +
	(energy >= 2.0) * (abs(eta) < 0.7)*(0.94) +
	(energy >= 2.0) * (abs(eta) >=0.7 && abs(eta) <=2.5) * (0.9)	}

}


##################
# Photon isolation
##################

module Isolation PhotonIsolation {
    set CandidateInputArray PhotonEfficiency/photons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray photons

    set DeltaRMax 0.1

    set PTMin 0.5

    set PTRatioMax 0.2
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
    set InputArray ElectronFilter/electrons
    set OutputArray electrons

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    set EfficiencyFormula {
	(energy < 3.0 ) * ( 0.00 ) +
    (abs(eta) > 2.50) * ( 0.00 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.58 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.7 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.6 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.7 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.8 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 0.69)                    * (0.84 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (  0.6 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.76 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.67 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.78 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.86 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 0.69)                    * ( 0.88 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * (  0.6 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (  0.8 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.68 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.84 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.88 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 0.69)                    * (  0.9 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.64 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.82 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.7 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.84 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.9 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 0.69)                    * (0.92 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * (0.64 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.86 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.74 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.87 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.91 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 0.69)                    * (0.94 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)   * (0.67 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.88 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.78 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.9 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.94 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 0.69)                    * (0.94 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.68 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.9 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.86 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.92 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.94 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 0.69)                    * (0.96 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (  0.7 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.92 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (  0.8 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.94 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.96 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 0.69)                    * ( 0.97 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * (0.68 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.96 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.84 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.94 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.98 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 0.69)                    * (0.98 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * ( 0.68 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.97 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.86 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.96 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.98 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 0.69)                    * ( 0.98 ) +
	( energy >=400  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.68 ) +
	( energy >=400  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.96 ) +
	( energy >=400  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.82 ) +
	( energy >=400  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.96 ) +
	( energy >=400  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.98 ) +
	( energy >=400  ) * (abs(eta) <= 0.69)                    * (0.98 )
    }
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
    set CandidateInputArray ElectronEfficiency/electrons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray electrons

    set DeltaRMax 0.1

    set PTMin 0.5

    set PTRatioMax 0.2
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
    set InputArray MuonFilter/muons
    set OutputArray muons

    # set EfficiencyFormula {efficiency as a function of eta and pt}
    # CONTROLLARE TAGLIO IN ENERGIA

    set EfficiencyFormula {
	(energy < 2.5 )     * (0.00) +
	(energy>=2.5  )     * (0.999)
    }
}

################
# Muon isolation
################

module Isolation MuonIsolation {
    set CandidateInputArray MuonEfficiency/muons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray muons

    set DeltaRMax 0.1

    set PTMin 0.5

    set PTRatioMax 0.2
}



###################
# Missing ET merger
###################

module Merger MissingET {
    # add InputArray InputArray
    add InputArray EFlowMerger/eflow
    set MomentumOutputArray momentum
}


##################
# Scalar HT merger
##################

module Merger ScalarHT {
    # add InputArray InputArray
    add InputArray EFlowMerger/eflow
    set EnergyOutputArray energy
}
######################
# EFlowFilter (UniqueObjectFinder)
######################
module UniqueObjectFinder EFlowFilter {
    add InputArray PhotonIsolation/photons photons
    add InputArray ElectronIsolation/electrons electrons
    add InputArray MuonIsolation/muons muons
    add InputArray EFlowMerger/eflow eflow
}

#################
# Neutrino Filter
#################

module PdgCodeFilter NeutrinoFilter {

    set InputArray Delphes/stableParticles
    set OutputArray filteredParticles

    set PTMin 0.0

    add PdgCode {12}
    add PdgCode {14}
    add PdgCode {16}
    add PdgCode {-12}
    add PdgCode {-14}
    add PdgCode {-16}

}

#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
    set InputArray NeutrinoFilter/filteredParticles

    set OutputArray jets

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 anti-kt with winner-take-all axis (for N-subjettiness), 8 N-jettiness, 9 Valencia
    set JetAlgorithm 9
    set ParameterR 0.5

    set JetPTMin 15.0
}

############
# Jet finder
############

module FastJetFinder FastJetFinderKt {
    #  set InputArray Calorimeter/towers
    set InputArray EFlowFilter/eflow

    set OutputArray KTjets

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 anti-kt with winner-take-all axis (for N-subjettiness), 8 N-jettiness, 9 Valencia
    set JetAlgorithm 4
    set ParameterR 0.5

    set JetPTMin 15.0
}

############
# Jet finder VLC
############
#R05 N2
# module FastJetFinder FastJetFinderVLC_R05_N2 {
#     #  set InputArray Calorimeter/towers
#     set InputArray EFlowFilter/eflow

#     set OutputArray VLCjetsR05N2

#     # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 anti-kt with winner-take-all axis (for N-subjettiness), 8 N-jettiness, 9 Valencia
#     set NJets 2
#     set ExclusiveClustering true
#     set JetAlgorithm 9
#     set ParameterR 0.5
#     set Beta 1.0
#     set Gamma 1.0

#     set JetPTMin 20.0
# }

########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray FastJetFinderKt/KTjets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

module BTagging BTaggingWP70 {
	set JetInputArray FastJetFinderKt/KTjets
	set BitNumber 0

	# efficiency formula for b-jets
    add EfficiencyFormula {5} {0.7}
    # default efficiency formula (misidentification rate)
    add EfficiencyFormula {0} {(energy >= 500 )*
        (energy >= 500 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
                    (energy >= 500 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 9e-3 )+ \
                    (energy >= 500 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 8e-3 )+ \
                    (energy >= 500 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 8e-3 )+ \
                    (energy >= 500 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
                    (energy >= 500 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 9e-3 )+ \
                    (energy >= 500 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
                    (energy >= 500 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-2 )+ \
                    (energy >= 500 )*( abs(eta) <=0.09                     ) * ( 1e-2 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 7e-3 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 6e-3 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 6e-3 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 5e-3 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 5e-3 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 5e-3 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 5e-3 )+ \
                    (energy < 500 && energy >= 250)*( abs(eta) <=0.09                     ) * ( 5e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 2e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 2e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 2e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 2e-3 )+ \
                    (energy < 250 && energy >= 100)*( abs(eta) <=0.09                     ) * ( 2e-3 )+ \
                    (energy < 100 )*( abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 2e-2 )+ \
                    (energy < 100 )*( abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 6e-3 )+ \
                    (energy < 100 )*( abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-3 )+ \
                    (energy < 100 )*( abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-3 )+ \
                    (energy < 100 )*( abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-3 )+ \
                    (energy < 100 )*( abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-3 )+ \
                    (energy < 100 )*( abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-3 )+ \
                    (energy < 100 )*( abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-3 )+ \
                    (energy < 100 )*( abs(eta) <=0.09                     ) * ( 1e-3 )
    }

    # efficiency formula for c-jets (misidentification rate)
    add EfficiencyFormula {4} {
        (energy >= 500 )*(abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
        (energy >= 500 )*(abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-2 )+ \
        (energy >= 500 )*(abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-2 )+ \
        (energy >= 500 )*(abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 3e-2 )+ \
        (energy >= 500 )*(abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 4e-2 )+ \
        (energy >= 500 )*(abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 4e-2 )+ \
        (energy >= 500 )*(abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 5e-2 )+ \
        (energy >= 500 )*(abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 5e-2 )+ \
        (energy >= 500 )*(abs(eta) <=0.09                     ) * ( 6e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 3e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 1e-2 )+ \
        (energy < 500 && energy >= 250)*(abs(eta) <=0.09                     ) * ( 1e-2 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 4e-2 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 2e-2 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 1e-2 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 8e-3 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 7e-3 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 7e-3 )+ \
        (energy < 250 && energy >= 100)*(abs(eta) <=0.09                     ) * ( 7e-3 )+ \
        (energy < 100 )* (abs(eta) <=2.66 && abs(eta) > 2.09  ) * ( 1e-1 )+ \
        (energy < 100 )* (abs(eta) <=2.09 && abs(eta) > 1.53  ) * ( 7e-2 )+ \
        (energy < 100 )* (abs(eta) <=1.53 && abs(eta)> 1.165  ) * ( 4e-2 )+ \
        (energy < 100 )* (abs(eta) <=1.165 && abs(eta)>0.885  ) * ( 2e-2 )+ \
        (energy < 100 )* (abs(eta) <=0.885 && abs(eta)>0.655  ) * ( 1e-2 )+ \
        (energy < 100 )* (abs(eta) <=0.655 && abs(eta)> 0.455 ) * ( 1e-2 )+ \
        (energy < 100 )* (abs(eta) <=0.455 && abs(eta)>0.27   ) * ( 1e-2 )+ \
        (energy < 100 )* (abs(eta) <=0.27 && abs(eta)> 0.09   ) * ( 9e-3 )+ \
        (energy < 100 )* (abs(eta) <=0.09                     ) * ( 9e-3 )
    }

}



##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass
    add Branch Delphes/allParticles Particle GenParticle

    ####

    # add Branch GenMissingET/momentum GenMissingET MissingET

    add Branch TrackMerger/tracks Track Track
    add Branch Calorimeter/towers Tower Tower

    add Branch HCal/eflowTracks EFlowTrack Track
    add Branch ECal/eflowPhotons EFlowPhoton Tower
    add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

    add Branch EFlowFilter/photons Photon Photon
    add Branch EFlowFilter/electrons Electron Electron
    add Branch EFlowFilter/muons Muon Muon

    add Branch GenJetFinder/jets GenJet Jet
    add Branch FastJetFinderKt/KTjets PFJet Jet

}
