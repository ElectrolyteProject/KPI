#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009,2010,2011,2012 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topolammps.tcl,v 1.44 2016/11/04 05:57:55 johns Exp $

# high level subroutines for LAMMPS support.
#
# import a LAMMPS data file.
# this behaves almost like a molfile plugin and will create a
# new molecule (but no representation) and return its molecule id.
#
# Arguments:
# filename = name of data file
# style = atomstyle
# flags = more flags. (currently not used)
proc readlammpsdata {filename style {flags none}} {
    global M_PI
    if {[catch {open $filename r} fp]} {
        vmdcon -err "readlammpsdata: problem opening data file: $fp\n"
        return -1
    }

    # parse lammps header section.
    array set lammps [readlammpsheader $fp]
    if {$lammps(atoms) < 1} {
        vmdcon -err "readlammpsdata: failed to parse lammps data header. abort."
        return -1
    }

    # create an empty molecule and timestep
    set mol -1
    if {[catch {mol new atoms $lammps(atoms)} mol]} {
        vmdcon -err "readlammpsdata: problem creating empty molecule: $mol"
        return -1
    } else {
        animate dup $mol
    }
    mol rename $mol [file tail $filename]
    set sel [atomselect $mol all]
    set boxdim {}

    if {($lammps(xy) != {}) && ($lammps(xz) != {}) && ($lammps(yz) != {})} {
        set lammps(triclinic) 1
        vmdcon -info "readlammpsdata: detected triclinic cell."
        set a [expr {$lammps(xhi) - $lammps(xlo)}]
        set ly [expr {$lammps(yhi) - $lammps(ylo)}]
        set lz [expr {$lammps(zhi) - $lammps(zlo)}]
        set b [expr {sqrt($ly*$ly + $lammps(xy)*$lammps(xy))}]
        set c [expr {sqrt($lz*$lz + $lammps(xz)*$lammps(xz)
                          + $lammps(yz)*$lammps(yz))}]
        set alpha [expr {($lammps(xy)*$lammps(xz) + $ly*$lammps(yz))/($b*$c)}]
        set beta  [expr {$lammps(xz)/$c}]
        set gamma [expr {$lammps(xy)/$b}]
        set alpha [expr {90.0 - asin($alpha)*180.0/$M_PI}]
        set beta  [expr {90.0 - asin($beta)*180.0/$M_PI}]
        set gamma [expr {90.0 - asin($gamma)*180.0/$M_PI}]

        set boxdim [list $a $b $c $alpha $beta $gamma]
        molinfo $mol set {a b c alpha beta gamma} $boxdim
        lappend boxdim $lammps(triclinic) $a $ly $lz $lammps(xy) $lammps(xz) $lammps(yz)
    } else {
        set $lammps(triclinic) 0
        set a [expr {$lammps(xhi) - $lammps(xlo)}]
        set b [expr {$lammps(yhi) - $lammps(ylo)}]
        set c [expr {$lammps(zhi) - $lammps(zlo)}]
        set boxdim [list $a $b $c 90.0 90.0 90.0]
        molinfo $mol set {a b c alpha beta gamma} $boxdim
        lappend boxdim $lammps(triclinic) $a $b $c 0.0 0.0 0.0
    }
    set atomidmap {}
    set atommasses {}

    # now loop through the file until we find a known section header.
    # then call a subroutine that parses this section. those subroutines
    # all take the current line number as last argument and return the
    # new value, so we can keep track of the current line and print more
    # useful error messages or warnings.
    set lineno $lammps(lineno)
    while {[gets $fp line] >= 0} {
        incr lineno
        if {[regexp {^\s*Atoms} $line ]} {
            # use atom style indicated by CGCMM header or use style hint comment.
            set stylehint $lammps(style)
            regexp {^\s*Atoms\s+#\s*([a-z]+)} $line x stylehint
            # for requested atom style 'auto' use the hint value instead
            if {[string equal $style auto]} { set style $stylehint }
            if {[string equal $style unknown]} {
                vmdcon -warn "readlammpsdata: automatic atom style detection requested,"
                vmdcon -warn "readlammpsdata: but no atom style hints in data file."
                vmdcon -warn "readlammpsdata: assuming atom style 'full' instead."
                set style {full}
            }
            # check for atom style consistency
            if {![string equal unknown $stylehint] && ![string equal $style $stylehint]} {
                vmdcon -warn "readlammpsdata: requested atom style '$style' is different from"
                vmdcon -warn "readlammpsdata: style hint '$stylehint' encoded into data file."
                vmdcon -warn "readlammpsdata: this may not work. Continuing with '$style'"
            }
            # if atom style is supported
            if { ![regexp {^(atomic|bond|angle|molecular|charge|full|sphere)} $style] } {
                vmdcon -err "readlammpsdata: unsupported atom style '$style'"
                return -1
            }
            set lineno [readlammpsatoms $fp $sel $style $lammps(cgcmm) $boxdim $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Atoms section."
                return -1
            }
            # retrieve map of atomids from user field and convert back to integer.
            set atomidmap {}
            foreach id [$sel get user] {
                lappend atomidmap [expr {int($id + 0.5)}]
            }
        } elseif {[regexp {^\s*Velocities} $line ]} {
            set lineno [readlammpsvelocities $fp $sel $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Velocities section."
                return -1
            }
        } elseif {[regexp {^\s*Masses} $line ]} {
            set lineno [readlammpsmasses $fp $mol $lammps(atomtypes) $lammps(cgcmm) atommasses $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Masses section."
                return -1
            }
        } elseif {[regexp {^\s*Bonds} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Bonds in data file"
                return -1
            }
            set lineno [readlammpsbonds $fp $sel $lammps(bonds) $atomidmap $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Bonds section."
                return -1
            }
        } elseif {[regexp {^\s*Angles} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Angles in data file"
                return -1
            }
            set lineno [readlammpsangles $fp $sel $lammps(angles) $atomidmap $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Angles section."
                return -1
            }
        } elseif {[regexp {^\s*Dihedrals} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Dihedrals in data file"
                return -1
            }
            set lineno [readlammpsdihedrals $fp $sel $lammps(dihedrals) $atomidmap $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Dihedrals section."
                return -1
            }
        } elseif {[regexp {^\s*Impropers} $line ]} {
            if {[llength $atomidmap] < 1} {
                vmdcon -err "readlammpsdata: Atoms section must come before Impropers in data file"
                return -1
            }
            set lineno [readlammpsimpropers $fp $sel $lammps(impropers) $atomidmap $lineno]
            if {$lineno < 0} {
                vmdcon -err "readlammpsdata: error reading Impropers section."
                return -1
            }
        } elseif {[regexp {^\s*(Pair Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(atomtypes) $lineno]
        } elseif {[regexp {^\s*(PairIJ Coeffs)} $line ]} {
            set skip [expr {$lammps(atomtypes)*($lammps(atomtypes)+1)/2}]
            set lineno [skiplammpslines $fp $skip $lineno]
        } elseif {[regexp {^\s*(Bond Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(bondtypes) $lineno]
        } elseif {[regexp {^\s*(Angle Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(angletypes) $lineno]
        } elseif {[regexp {^\s*(BondBond Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(angletypes) $lineno]
        } elseif {[regexp {^\s*(BondAngle Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(angletypes) $lineno]
        } elseif {[regexp {^\s*(Dihedral Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(MiddleBondTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(EndBondTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(AngleTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(AngleAngleTorsion Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(BondBond13 Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(dihedraltypes) $lineno]
        } elseif {[regexp {^\s*(Improper Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(impropertypes) $lineno]
        } elseif {[regexp {^\s*(AngleAngle Coeffs)} $line ]} {
            set lineno [skiplammpslines $fp $lammps(impropertypes) $lineno]
        } elseif { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines silently
        } else {
            vmdcon -err "readlammpsdata: unkown header line: $lineno : $line "
            vmdcon -err "readlammpsdata: cannot continue. "
            return -1
        }
        set lammps(lineno) $lineno
    }
    close $fp

    # apply masses. Atoms section sets a default of 1.0.
    # since the Masses section can appear before the Atoms section
    # we have to set it here after the parsing.
    if {[llength atommasses] > 0} {
        foreach {t m} $atommasses {
            set msel [atomselect $mol "type '$t'"]
            $msel set mass $m
            $msel delete
        }
    }

    mol reanalyze $mol
    #variable newaddsrep
    #if {$newaddsrep} {
    #    adddefaultrep $mol
    #}
    #$sel delete

    return $mol
}

# skip over a given number of non-empty, non-comment lines
proc skiplammpslines {fp num lineno} {

    while {[gets $fp line] >= 0} {
        if {$num <= 0} break
        incr lineno
        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
          incr num -1
        }
    }
    return $lineno
}

# read lammps header section from opened file
# and return as an array.
proc readlammpsheader {fp} {
    array set lammps {
        atoms 0 atomtypes 0 bonds 0 bondtypes 0 angles 0 angletypes 0
        dihedrals 0 dihedraltypes 0 impropers 0 impropertypes 0 xtrabond 0
        xlo -0.5 xhi 0.5 ylo -0.5 yhi 0.5 zlo -0.5 zhi 0.5 xy {} xz {} yz {}
        lineno 0 cgcmm 0 triclinic 0 style unknown typemass 1
    }
    set x {}

    vmdcon -info "parsing LAMMPS header."

    # first header line is skipped by LAMMPS. so we put a flag to
    # detect CMM style CG data files with additional information.
    gets $fp line
    if {[string match "*CGCMM*" $line]} {
        set lammps(cgcmm) 1
        vmdcon -info "detected CGCMM style file. will try to parse additional data."
        if {[string match "*atom_style*" $line]} {
            if { [regexp {^.*atom_style\s+(atomic|bond|angle|molecular|charge|full|sphere)\s*.*}  $line x lammps(style) ] } {
                vmdcon -info "Probable atom_style: $lammps(style)"
            }
        }
    }
    set lineno 1
    set offs [tell $fp]
    set lammps(lineno) $lineno

    #
    while {[gets $fp line] >= 0} {
        incr lineno
        if { [      regexp {^\s*(\d+)\s+atoms}     $line x lammps(atoms) ] } {
        } elseif { [regexp {^\s*(\d+)\s+bonds}     $line x lammps(bonds) ] } {
        } elseif { [regexp {^\s*(\d+)\s+angles}    $line x lammps(angles)] } {
        } elseif { [regexp {^\s*(\d+)\s+dihedrals} $line x lammps(dihedrals)] } {
        } elseif { [regexp {^\s*(\d+)\s+impropers} $line x lammps(impropers)] } {
        } elseif { [regexp {^\s*(\d+)\s+atom types}     $line x lammps(atomtypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+bond types}     $line x lammps(bondtypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+angle types}    $line x lammps(angletypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+dihedral types} $line x lammps(dihedraltypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+improper types} $line x lammps(impropertypes)] } {
        } elseif { [regexp {^\s*(\d+)\s+extra bond per atom} $line x lammps(xtrabond)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+xlo xhi} $line \
                        x lammps(xlo) lammps(xhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+ylo yhi} $line \
                        x lammps(ylo) lammps(yhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+zlo zhi} $line \
                        x lammps(zlo) lammps(zhi)] } {
        } elseif { [regexp {^\s*([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+([-[:digit:].Ee+]+)\s+xy\s+xz\s+yz} $line x lammps(xy) lammps(xz) lammps(yz)] } {
        } elseif { [regexp {^\s*(\#.*|)$} $line ] } {
        } elseif {[regexp {^\s*(Atoms|Velocities|Masses|Shapes|Dipoles|Bonds|Angles|Dihedrals|Impropers|(Pair|Bond|Angle|Dihedral|Improper) Coeffs)} $line ]} {
            seek $fp $offs start
            break
        } else {
            vmdcon -warn "readlammpsheader: skipping unkown header line: $lineno : $line "
        }
        set offs [tell $fp]
        set lammps(lineno) $lineno
    }

    return [array get lammps]
}

# parse atom section
proc readlammpsatoms {fp sel style cgcmm boxdata lineno} {
    global M_PI
    set numatoms [$sel num]
    set atomdata {}
    set boxx 0.0
    set boxy 0.0
    set boxz 0.0
    set alpha 90.0
    set beta  90.0
    set gamma 90.0
    set triclinic 0
    set lx 0.0
    set ly 0.0
    set lz 0.0
    set xy 0.0
    set xz 0.0
    set yz 0.0

    lassign $boxdata boxx boxy boxz alpha beta gamma triclinic lx ly lz xy xz yz

    vmdcon -info "parsing LAMMPS Atoms section with style '$style'."

    set curatoms 0
    while {[gets $fp line] >= 0} {
        incr lineno

        set atomid 0
        set resid 0
        set atomname ""
        set resname ""
        set atomtype 0
        set charge 0.0
        set mass 1.0 ; #  m=0.0 in MD gets us in trouble, so use a different default.
        set radius 1.5 ; # default radius for unknown elements.
        # XXX: we could have a guess(element|radius) utility for setting this to something better.
        set x 0.0
        set y 0.0
        set z 0.0
        set xi 0
        set yi 0
        set zi 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty, whitespace or comment lines.
        } else {
            if {[regexp {^(.*)\#\s*(\S+)(\s+(\S+))?} $line all nline atomname dummy resname]} {
                set line $nline
            }
            incr curatoms
            switch $style { # XXX: use regexp based parser to detect wrong formats.

                atomic {
                    if {[llength $line] >= 8} {
                        lassign $line atomid       atomtype        x y z xi yi zi
                    } else {
                        lassign $line atomid       atomtype        x y z
                    }
                }

                bond      -
                angle     -
                molecular {
                    if {[llength $line] >= 9} {
                        lassign $line atomid resid atomtype        x y z xi yi zi
                    } else {
                        lassign $line atomid resid atomtype        x y z
                    }
                }

                charge {
                    if {[llength $line] >= 9} {
                        lassign $line atomid atomtype charge x y z xi yi zi
                    } else {
                        lassign $line atomid atomtype charge x y z
                    }
                }

                sphere {
                    if {[llength $line] >= 10} {
                        lassign $line atomid atomtype radius mass x y z xi yi zi
                    } else {
                        lassign $line atomid atomtype radius mass x y z
                    }
                    # sphere has diameter and density instead of radius and mass
                    # convert them accordingly
                    set radius [expr {0.5*$radius}]
                    set mass [expr {4.0/3.0*$M_PI*$radius*$radius*$radius*$mass}]
                }

                full {
                    if {[llength $line] >= 10} {
                        lassign $line atomid resid atomtype charge x y z xi yi zi
                    } else {
                        lassign $line atomid resid atomtype charge x y z
                    }
                }

                default   {
                    # ignore this unsupported style
                }
            }

            # sanity check on input data. if x/y/z are empty the atom style
            # must have been specified wrong. sufficient to check for $z
            if { [string length $z] == 0 } {
                vmdcon -err "readlammpsatoms: not enough data for style '$style' in line $lineno."
                return -1
            }

            # XXX: no reason to discriminate by atomid. we store the
            # atomid in the user field and assign atoms sorted by atomid
            # (for faster lookup) and can retrieve that info later and
            # use it for mapping bonds, angles, etc to it.
            ################
            # if {$atomid > $numatoms} {
            #     vmdcon -err "readlammpsatoms: only atomids 1-$numatoms are supported. $lineno : $line "
            #    return -1
            # }
            if {$cgcmm} {
                if {[string length $atomname]} {
                    set atomtype $atomname ; # if we have CGCMM data use that.
                } else {
                    set atomname $atomtype
                }
            } else {
                set atomname $atomtype
            }
            if {$triclinic} {
                lappend atomdata [list $atomid $resid $resname $atomname $atomtype $charge \
                                      [expr {$x + $xi*$lx + $yi*$xy + $zi*$xz}] \
                                      [expr {$y + $yi*$ly + $zi*$yz}] \
                                      [expr {$z + $zi*$lz}] $mass $radius ]
            } else {
                lappend atomdata [list $atomid $resid $resname $atomname $atomtype $charge \
                                      [expr {$xi*$boxx + $x}] [expr {$yi*$boxy + $y}] \
                                      [expr {$zi*$boxz + $z}] $mass $radius ]
            }
        }
        if {$curatoms >= $numatoms} break
    }
    vmdcon -info "applying atoms data. sorted by atom id."
    $sel set {user resid resname name type charge x y z mass radius} \
        [lsort -integer -index 0 $atomdata]
    return $lineno
}

# parse masses section
proc readlammpsmasses {fp mol numtypes cgcmm massmap lineno} {
    vmdcon -info "parsing LAMMPS Masses section."

    upvar $massmap massdata
    set massdata {}
    set curtypes 0
    while {[gets $fp line] >= 0} {
        incr lineno

        set typeid 0
        set mass 0.0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curtypes
            set typename {}
            if {[regexp {^(.*)\#\s*(\S+).*} $line all nline typename]} {
                set line $nline
            }
            lassign $line typeid mass
            if {$typeid > $numtypes} {
                vmdcon -err "readlammpsmasses: only typeids 1-$numtypes are supported. $lineno : $line "
                return -1
            }
            # if we have a CGCMM style data file, we have strings for types.
            if {$cgcmm && ([string length $typename] > 0)} {
                lappend massdata $typename $mass
            } else {
                lappend massdata $typeid $mass
            }
        }
        if {$curtypes >= $numtypes} break
    }
    return $lineno
}

# parse velocities section
proc readlammpsvelocities {fp sel lineno} {
    set numatoms [$sel num]
    set velocitydata {}

    vmdcon -info "parsing LAMMPS Velocities section."

    set curatoms 0
    set veldata {}
    while {[gets $fp line] >= 0} {
        incr lineno

        set atomid 0
        set vx 0.0
        set vy 0.0
        set vz 0.0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curatoms
            lassign $line atomid vx vy vz

            #if {$atomid > $numatoms} {
            #    vmdcon -err "readlammpsvelocities: only atomids 1-$numatoms are supported. $lineno : $line "
            #    return -1
            #}
            lappend veldata [list $atomid $vx $vy $vz]
        }
        if {$curatoms >= $numatoms} break
    }
    if { [catch {$sel set {user vx vy vz} [lsort -integer -index 0 $veldata]} errmsg] } {
        vmdcon -warn "readlammpsvelocities: problems assigning velocities. skipping..."
    }
    return $lineno
}


# parse bond section
proc readlammpsbonds {fp sel numbonds atomidmap lineno} {
    set curbonds 0
    set bonddata {}

    vmdcon -info "parsing LAMMPS Bonds section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curbonds
            lassign $line num type a b ;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            if {($aidx < 0) || ($bidx < 0)} {
                vmdcon -err "readlammpsdata: bond with non-existent atomid on line: $lineno"
                return -1
            }
            lappend bonddata [list $aidx $bidx $type]
        }
        if {$curbonds >= $numbonds} break
    }
    vmdcon -info "applying bonds data."
    setbondlist $sel type $bonddata
    return $lineno
}

# parse angle section
proc readlammpsangles {fp sel numangles atomidmap lineno} {
    set curangles 0
    set angledata {}

    vmdcon -info "parsing LAMMPS Angles section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set c 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curangles
            lassign $line num type a b c ;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            set cidx [lsearch -sorted -integer $atomidmap $c]
            if {($aidx < 0) || ($bidx < 0) || ($cidx < 0)} {
                vmdcon -err "readlammpsdata: angle with non-existent atomid on line: $lineno"
                return -1
            }
            lappend angledata [list $type $aidx $bidx $cidx]
        }
        if {$curangles >= $numangles} break
    }
    vmdcon -info "applying angles data."
    setanglelist $sel $angledata
    return $lineno
}

# parse dihedral section
proc readlammpsdihedrals {fp sel numdihedrals atomidmap lineno} {
    set curdihedrals 0
    set dihedraldata {}

    vmdcon -info "parsing LAMMPS Dihedrals section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set c 0
        set d 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curdihedrals
            lassign $line num type a b c d ;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            set cidx [lsearch -sorted -integer $atomidmap $c]
            set didx [lsearch -sorted -integer $atomidmap $d]
            if {($aidx < 0) || ($bidx < 0) || ($cidx < 0) || ($didx < 0)} {
                vmdcon -err "readlammpsdata: dihedral with non-existent atomid on line: $lineno"
                return -1
            }
            lappend dihedraldata [list $type $aidx $bidx $cidx $didx]
        }
        if {$curdihedrals >= $numdihedrals} break
    }
    vmdcon -info "applying dihedrals data."
    setdihedrallist $sel $dihedraldata
    return $lineno
}

# parse improper section
proc readlammpsimpropers {fp sel numimpropers atomidmap lineno} {
    set curimpropers 0
    set improperdata {}

    vmdcon -info "parsing LAMMPS Impropers section."
    while {[gets $fp line] >= 0} {
        incr lineno

        set num 0
        set a 0
        set b 0
        set c 0
        set d 0
        set type 0

        if { [regexp {^\s*(\#.*|)$} $line ] } {
            # skip empty lines.
        } else {
            incr curimpropers
            lassign $line num type a b c d;# XXX: use regexp based parser to detect wrong format.
            # map atomid to vmd atom index.
            set aidx [lsearch -sorted -integer $atomidmap $a]
            set bidx [lsearch -sorted -integer $atomidmap $b]
            set cidx [lsearch -sorted -integer $atomidmap $c]
            set didx [lsearch -sorted -integer $atomidmap $d]
            if {($aidx < 0) || ($bidx < 0) || ($cidx < 0) || ($didx < 0)} {
                vmdcon -err "readlammpsdata: improper with non-existent atomid on line: $lineno"
                return -1
            }
            lappend improperdata [list $type $aidx $bidx $cidx $didx]
        }
        if {$curimpropers >= $numimpropers} break
    }
    vmdcon -info "applying impropers data."
    setimproperlist $sel $improperdata
    return $lineno
}


# export internal structure data to a LAMMPS data file.
# this requires that a corresponding set of information
# is already present in VMD's memory.
# Arguments:
# mol = molecule id with matching coordinate data
# filename = name of data file
# style = atom style
# sel = selection function for the subset to be written out.
# flags = more flags. (currently not used)
proc writelammpsdata {mol filename style sel {flags none}} {
    if {[catch {open $filename w} fp]} {
        vmdcon -err "writelammpsdata: problem opening data file: $fp\n"
        return -1
    }

    # initialize system default settings
    array set lammps {
        atoms 0 atomtypes 0 bonds 0 bondtypes 0 angles 0 angletypes 0
        dihedrals 0 dihedraltypes 0 impropers 0 impropertypes 0 xtrabond 0
        xlo 0 xhi 0 ylo 0 yhi 0 zlo 0 zhi 0 xy 0 xz 0 yz 0 triclinic 0
        style unknown typemass 1
    }

    # gather available system information
    set lammps(atoms)         [$sel num]
    set lammps(bonds)         [bondinfo     numbonds     $sel]
    set lammps(angles)        [angleinfo    numangles    $sel]
    set lammps(dihedrals)     [dihedralinfo numdihedrals $sel]
    set lammps(impropers)     [improperinfo numimpropers $sel]
    set lammps(atomtypes)     [llength [lsort -ascii -unique [$sel get {type}]]]
    set lammps(bondtypes)     [bondinfo     numbondtypes     $sel]
    set lammps(angletypes)    [angleinfo    numangletypes    $sel]
    set lammps(dihedraltypes) [dihedralinfo numdihedraltypes $sel]
    set lammps(impropertypes) [improperinfo numimpropertypes $sel]
    set lammps(style) $style

    # correct system information to allow only information valid
    # for the selected atom style
    switch $style {
        atomic -
        charge {
            set lammps(bonds) 0
            set lammps(angles) 0
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(bondtypes) 0
            set lammps(angletypes) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        sphere {
            set lammps(typemass) 0
            set lammps(bonds) 0
            set lammps(angles) 0
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(bondtypes) 0
            set lammps(angletypes) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        bond  {
            set lammps(angles) 0
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(angletypes) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        angle  {
            set lammps(dihedrals) 0
            set lammps(impropers) 0
            set lammps(dihedraltypes) 0
            set lammps(impropertypes) 0
        }

        molecular -
        full -
        default { ; # print all sections
        }
    }

    # initialize simulation cell dimensions from min/max search
    lassign [measure minmax $sel -withradii] min max
    lassign $min xlo ylo zlo
    lassign $max xhi yhi zhi
    lassign [molinfo $mol get {a b c alpha beta gamma}] boxx boxy boxz alpha beta gamma

    # override min/max settings where box info available.
    # try to (mostly) preserve the center of the cell, by
    # deriving it from the preceeding min/max search.
    set small 0.0001
    if {$boxx > $small} {
        set lammps(xmid) [expr {($xlo + $xhi) * 0.5}]
        set lammps(xlo)  [expr {-0.5*$boxx + $lammps(xmid)}]
        set lammps(xhi)  [expr { 0.5*$boxx + $lammps(xmid)}]
    } else {
        set lammps(xmid)  0.0
        set lammps(xlo)  -0.5
        set lammps(xhi)   0.5
    }
    if {$boxy > $small} {
        set lammps(ymid) [expr {($ylo + $yhi) * 0.5}]
        set lammps(ylo)  [expr {-0.5*$boxy + $lammps(ymid)}]
        set lammps(yhi)  [expr { 0.5*$boxy + $lammps(ymid)}]
    } else {
        set lammps(ymid)  0.0
        set lammps(ylo)  -0.5
        set lammps(yhi)   0.5
    }
    if {$boxz > $small} {
        set lammps(zmid) [expr {($zlo + $zhi) * 0.5}]
        set lammps(zlo)  [expr {-0.5*$boxz + $lammps(zmid)}]
        set lammps(zhi)  [expr { 0.5*$boxz + $lammps(zmid)}]
    } else {
        set lammps(zmid)  0.0
        set lammps(zlo)  -0.5
        set lammps(zhi)   0.5
    }

    # if angle is not set assume orthogonal.
    if {$alpha == 0.0} { set alpha 90.0 }
    if {$beta  == 0.0} { set beta  90.0 }
    if {$gamma == 0.0} { set gamma 90.0 }

    if {($alpha != 90.0) || ($beta != 90.0) || ($gamma != 90.0)} {
        set conv 0.01745329251994329576; # pi/180.0
        set lammps(triclinic) 1
        set lammps(xy) [expr {$boxy * cos($gamma*$conv)}]
        set lammps(xz) [expr {$boxz * cos($beta*$conv)}]
        set lammps(yz) 0.0
        set ly [expr {sqrt($boxy*$boxy - $lammps(xy)*$lammps(xy))}]
        if {abs($ly) > $small} {
            set lammps(yz) [expr {($boxy*$boxz*cos($alpha*$conv)
                                   - $lammps(xy)*$lammps(xz)) / $ly}]
        }
        set lz [expr {sqrt($boxz*$boxz - $lammps(xz)*$lammps(xz) - $lammps(yz)*$lammps(yz))}]
        # update y/z-box boundaries for tilt
        if {$ly > $small} {
            set lammps(ylo)  [expr {-0.5*$ly + $lammps(ymid)}]
            set lammps(yhi)  [expr { 0.5*$ly + $lammps(ymid)}]
        }
        if {$lz > $small} {
            set lammps(zlo)  [expr {-0.5*$lz + $lammps(zmid)}]
            set lammps(zhi)  [expr { 0.5*$lz + $lammps(zmid)}]
        }
    }

    # write out supported data file sections
    writelammpsheader $fp [array get lammps]

    # write out hints about type to number mappings
    # for coefficient settings
    writelammpscoeffhint $fp $sel atoms
    if {$lammps(bonds) > 0} {
        writelammpscoeffhint $fp $sel bonds
    }
    if {$lammps(angles) > 0} {
        writelammpscoeffhint $fp $sel angles
    }
    if {$lammps(dihedrals) > 0} {
        writelammpscoeffhint $fp $sel dihedrals
    }
    if {$lammps(impropers) > 0} {
        writelammpscoeffhint $fp $sel impropers
    }
    if {$lammps(typemass) > 0} {
        writelammpsmasses $fp $sel
    }
    writelammpsatoms $fp $sel $style
    set atomidmap  [$sel list]
    if {$lammps(bonds) > 0} {
        writelammpsbonds $fp $sel $atomidmap
    }
    if {$lammps(angles) > 0} {
        writelammpsangles $fp $sel $atomidmap
    }
    if {$lammps(dihedrals) > 0} {
        writelammpsdihedrals $fp $sel $atomidmap
    }
    if {$lammps(impropers) > 0} {
        writelammpsimpropers $fp $sel $atomidmap
    }
    close $fp
    return 0
}


# write lammps header section to open file
proc writelammpsheader {fp flags} {
    variable version
    array set lammps $flags
    # first header line is skipped.
    puts $fp "LAMMPS data file. CGCMM style. atom_style $lammps(style) generated by VMD/TopoTools v$version on [clock format [clock seconds]]"

    foreach key {atoms bonds angles dihedrals impropers} {
        puts $fp [format " %d %s" $lammps($key) $key]
    }

    foreach key {atomtypes bondtypes  angletypes dihedraltypes impropertypes} {
        puts $fp [format " %d %s" $lammps($key) [regsub types $key " &"]]
    }

    puts $fp [format " %.6f %.6f  xlo xhi" $lammps(xlo) $lammps(xhi)]
    puts $fp [format " %.6f %.6f  ylo yhi" $lammps(ylo) $lammps(yhi)]
    puts $fp [format " %.6f %.6f  zlo zhi" $lammps(zlo) $lammps(zhi)]

    if {$lammps(triclinic)} {
        puts $fp [format " %.6f %.6f %.6f xy xz yz" $lammps(xy) $lammps(xz) $lammps(yz)]
    }

    puts $fp ""
    return
}

# write masses section, but only if number of masses
# matches the number of atom types and if no mass is < 0.01
proc writelammpsmasses {fp sel} {

    # first run the checks and build list of masses
    set typemap  [lsort -unique -ascii [$sel get type]]
    set masslist {}
    set mol [$sel molid]
    set selstr [$sel text]
    foreach type $typemap {
        set tsel [atomselect $mol "( $selstr ) and (type '$type')"]
        set mass [lsort -unique -real [$tsel get mass]]
        $tsel delete
        if {[llength $mass] != 1} return
        if {$mass < 0.01} return
        lappend masslist $mass
    }

    # we passed the test, write out what we learned.
    vmdcon -info "writing LAMMPS Masses section."

    puts $fp " Masses\n"
    set typeid 1
    foreach mass $masslist type $typemap {
        puts $fp [format " %d %.6f \# %s" $typeid $mass $type]
        incr typeid
    }
    puts $fp ""
    return
}

# write atoms section
proc writelammpsatoms {fp sel style} {
    global M_PI

    vmdcon -info "writing LAMMPS Atoms section in style '$style'."

    puts $fp " Atoms # $style\n"
    set typemap [lsort -unique -ascii [$sel get type]]
    set resmap  [lsort -unique -integer [$sel get resid]]
    set atomid 0
    foreach adat [$sel get {type resid charge x y z resname mass radius}] {
        lassign $adat type resid charge x y z resname mass radius
        set atomtype [lsearch -sorted -ascii $typemap $type]
        set resid    [lsearch -sorted -integer $resmap $resid]
        incr atomid
        incr atomtype
        incr resid
        switch $style {
            atomic {
                puts $fp [format "%d %d %.6f %.6f %.6f \# %s" \
                              $atomid        $atomtype  $x $y $z $type]
            }

            bond  -
            angle -
            molecular {
                puts $fp [format "%d %d %d %.6f %.6f %.6f \# %s %s" \
                              $atomid $resid $atomtype  $x $y $z $type $resname]
            }

            charge    {
                puts $fp [format "%d %d %.6f %.6f %.6f %.6f \# %s" \
                              $atomid $atomtype $charge $x $y $z $type]
            }

            sphere {
                # sphere has diameter and density instead of radius and mass
                # convert them accordingly
                set mass [expr {$mass/(4.0/3.0*$M_PI*$radius*$radius*$radius)}]
                set radius [expr {2.0*$radius}]
                puts $fp [format "%d %d %.6f %.6f %.6f %.6f %.6f \# %s" \
                              $atomid $atomtype $radius $mass $x $y $z $type]
            }

            full      {
                puts $fp [format "%d %d %d %.6f %.6f %.6f %.6f \# %s %s" \
                              $atomid $resid $atomtype $charge $x $y $z $type $resname]
            }

            default   {
                # ignore this unsupported style
                # XXX: add a way to flag an error. actually the test for
                #      supported lammps atom styles should be done on a
                #      much higher level, so that we don't do unneeded work.
            }
        }
    }
    puts $fp ""
    return
}

# write bond section
proc writelammpsbonds {fp sel atomidmap} {
    set bonddata  [bondinfo getbondlist   $sel type]
    set bondtypes [bondinfo bondtypenames $sel type]
    vmdcon -info "writing LAMMPS Bonds section."
    puts $fp " Bonds\n"

    set bondid 0
    foreach bdat $bonddata {
        incr bondid
        lassign $bdat a b t
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        set type [lsearch -ascii $bondtypes $t]

        # go from 0-based to 1-based indexing and write out
        incr at1; incr at2; incr type
        puts $fp [format "%d %d %d %d" $bondid $type $at1 $at2]
    }
    puts $fp ""
    return
}

# write angle section
proc writelammpsangles {fp sel atomidmap} {
    set angledata  [angleinfo getanglelist   $sel]
    set angletypes [angleinfo angletypenames $sel]
    vmdcon -info "writing LAMMPS Angles section."
    puts $fp " Angles\n"

    set angleid 0
    foreach adat $angledata {
        incr angleid
        lassign $adat t a b c
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        set at3  [lsearch -integer -sorted $atomidmap $c]
        set type [lsearch -ascii $angletypes $t]

        # go from 0-based to 1-based indexing and write out
        incr at1; incr at2; incr at3; incr type
        puts $fp [format "%d %d %d %d %d" $angleid $type $at1 $at2 $at3]
    }
    puts $fp ""
    return
}

# write dihedral section
proc writelammpsdihedrals {fp sel atomidmap} {
    set dihedraldata  [dihedralinfo getdihedrallist   $sel]
    set dihedraltypes [dihedralinfo dihedraltypenames $sel]
    vmdcon -info "writing LAMMPS Dihedrals section."
    puts $fp " Dihedrals\n"

    set dihedralid 0
    foreach adat $dihedraldata {
        incr dihedralid
        lassign $adat t a b c d
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        set at3  [lsearch -integer -sorted $atomidmap $c]
        set at4  [lsearch -integer -sorted $atomidmap $d]
        set type [lsearch -ascii $dihedraltypes $t]

        # go from 0-based to 1-based indexing and write out
        incr at1; incr at2; incr at3; incr at4; incr type
        puts $fp [format "%d %d %d %d %d %d" $dihedralid $type $at1 $at2 $at3 $at4]
    }
    puts $fp ""
    return
}

# write improper section
proc writelammpsimpropers {fp sel atomidmap} {
    set improperdata  [improperinfo getimproperlist   $sel]
    set impropertypes [improperinfo impropertypenames $sel]
    vmdcon -info "writing LAMMPS Impropers section."
    puts $fp " Impropers\n"

    set improperid 0
    foreach adat $improperdata {
        incr improperid
        lassign $adat t a b c d
        set at1  [lsearch -integer -sorted $atomidmap $a]
        set at2  [lsearch -integer -sorted $atomidmap $b]
        set at3  [lsearch -integer -sorted $atomidmap $c]
        set at4  [lsearch -integer -sorted $atomidmap $d]
        set type [lsearch -ascii $impropertypes $t]

        # go from 0-based to 1-based indexing and write out
        incr at1; incr at2; incr at3; incr at4; incr type
        puts $fp [format "%d %d %d %d %d %d" $improperid $type $at1 $at2 $at3 $at4]
    }
    puts $fp ""
    return
}

# returns 0 if lammps atom style is supported by topotools and 1 if not.
proc checklammpsstyle {style} {
    switch $style {

        auto -
        atomic -
        bond  -
        angle -
        molecular -
        charge -
        sphere -
        full {
            return 0
        }

        default {
            return 1
        }
    }
}

# write hints about type coefficient mappings
proc writelammpscoeffhint {fp sel type} {
    switch $type {
        atoms {
            puts $fp "\# Pair Coeffs\n\#"
            set aid 1
            set atlist [lsort -ascii -unique [$sel get {type}]]
            foreach at $atlist {
                puts $fp "\# $aid  $at"
                incr aid
            }
        }
        bonds {
            puts $fp "\# Bond Coeffs\n\#"
            set bid 1
            foreach bt [bondinfo bondtypenames $sel type] {
                puts $fp "\# $bid  $bt"
                incr bid
            }
        }
        angles {
            puts $fp "\# Angle Coeffs\n\#"
            set aid 1
            foreach at [angleinfo angletypenames $sel] {
                puts $fp "\# $aid  $at"
                incr aid
            }
        }
        dihedrals {
            puts $fp "\# Dihedral Coeffs\n\#"
            set did 1
            foreach dt [dihedralinfo dihedraltypenames $sel] {
                puts $fp "\# $did  $dt"
                incr did
            }
        }
        impropers {
            puts $fp "\# Improper Coeffs\n\#"
            set iid 1
            foreach it [improperinfo impropertypenames $sel] {
                puts $fp "\# $iid  $it"
                incr iid
            }
        }
        default {
            vmdcon -warn "writelammpscoeffhint: don't know how to write hints for '$type'"
            return 1
        }
    }
    puts $fp ""
    return
}

#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009,2010,2011 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topobonds.tcl,v 1.15 2014/08/19 16:45:04 johns Exp $

# Return info about bonds.
# we list and count only bonds that are entirely within the selection.
proc bondinfo {infotype sel {flag none}} {

    set numbonds 0
    set bidxlist {}
    array set bondtypes {}

    set aidxlist [$sel list]
    set bondlist [$sel getbonds]
    set btyplist [$sel getbondtypes]
    set bordlist [$sel getbondorders]

    foreach a $aidxlist bl $bondlist tl $btyplist ol $bordlist {
        foreach b $bl t $tl o $ol {
            if {($a < $b) && ([lsearch -sorted -integer $aidxlist $b] != -1)} {
                incr numbonds
                switch $flag {
                    type   {lappend bidxlist [list $a $b $t]}
                    order  {lappend bidxlist [list $a $b $o]}
                    both   {lappend bidxlist [list $a $b $t $o]}
                    lammps {lappend bidxlist [list $numbonds $a $b $t]}
                    none   {lappend bidxlist [list $a $b]}
                }
            }
            set bondtypes($t) 1
        }
    }

    switch $infotype {
        numbonds      { return $numbonds }
        numbondtypes  { return [array size bondtypes] }
        bondtypenames { return [lsort -ascii [array names bondtypes]] }
        getbondlist   { return $bidxlist }
        default       { return "bug? shoot the programmer!"}
    }
}

# delete all contained bonds of the selection.
proc clearbonds {sel} {

    # special optimization for "all" selection.
    if {[string equal "all" [$sel text]]} {
        set nulllist {}
        for {set i 0} {$i < [$sel num]} {incr i} {
            lappend nullist {}
        }
        $sel setbonds $nullist
        return
    }

    set mol [$sel molid]
    foreach b [bondinfo getbondlist $sel none] {
        delbond $mol [lindex $b 0] [lindex $b 1]
    }
}

# guess bonds from atom radii. Interface to "mol bondsrecalc".
# XXX: currently only works for selection "all".
proc guessbonds {sel} {

    set mol [$sel molid]
    # special optimization for "all" selection.
    if {[string equal "all" [$sel text]]} {
        # Use VMD's built-in bond determination heuristic to guess the bonds
        mol bondsrecalc $mol

        # Mark the bonds as "validated" so VMD will write
        # them out when the structure gets written out,
        # e.g. to a PSF file, even if no other bond editing was done.
        mol dataflag $mol set bonds

        return
    } else {
        vmdcon -err "topo guessbonds: this feature currently only works with an 'all' selection"
        return
    }
}

# reset bonds to data in bondlist
proc setbondlist {sel flag bondlist} {

    clearbonds $sel
    set nbnd [llength $bondlist]
    if {$nbnd == 0} { return 0}
    # set defaults
    set n 0
    set t unknown
    set o 1
    set mol [$sel molid]
    set a -1
    set b -1
    set fract  [expr {100.0/$nbnd}]
    set deltat 2000
    set newt   $deltat

    # special optimization for "all" selection.
    if {[string equal "all" [$sel text]]} {
        set nulllist {}
        for {set i 0} {$i < [$sel num]} {incr i} {
            set blist($i) $nulllist
            set olist($i) $nulllist
            set tlist($i) $nulllist
        }
        foreach bond $bondlist {
            switch $flag {
                type   {lassign $bond a b t  }
                order  {lassign $bond a b o  }
                both   {lassign $bond a b t o}
                lammps {lassign $bond n a b t}
                none   {lassign $bond a b    }
            }
            lappend blist($a) $b
            lappend blist($b) $a
            lappend olist($a) $o
            lappend olist($b) $o
            lappend tlist($a) $t
            lappend tlist($b) $t
        }
        set dlist {}
        for {set i 0} {$i < [$sel num]} {incr i} {
            lappend dlist $blist($i)
        }
        $sel setbonds $dlist
        set dlist {}
        for {set i 0} {$i < [$sel num]} {incr i} {
            lappend dlist $olist($i)
        }
        $sel setbondorders $dlist
        set dlist {}
        for {set i 0} {$i < [$sel num]} {incr i} {
            lappend dlist $tlist($i)
        }
        $sel setbondtypes $dlist
        return 0
    }

    # XXX: fixme!
    # using addbond is very inefficient with a large number of bonds
    # that are being added. it is better to fill the corresponding
    # bondlists directly. the code above should be better, but uses
    # much more memory and needs to be generalized.

    # XXX: add sanity check on data format
    set i 0
    foreach bond $bondlist {
        incr i
        set time [clock clicks -milliseconds]
        if {$time > $newt} {
            set percent [format "%3.1f" [expr {$i*$fract}]]
            vmdcon -info "setbondlist: $percent% done."
            display update ui
            set newt [expr {$time + $deltat}]
        }
        switch $flag {
            type   {lassign $bond a b t  }
            order  {lassign $bond a b o  }
            both   {lassign $bond a b t o}
            lammps {lassign $bond n a b t}
            none   {lassign $bond a b    }
        }
        addbond $mol $a $b $t $o
    }
    return 0
}

# guess bonds type names from atom types.
proc retypebonds {sel} {

    set bondlist  [bondinfo getbondlist $sel none]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newbonds {}

    foreach bond $bondlist {
        set idx [lsearch -sorted -integer $atomindex [lindex $bond 0]]
        set a [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex [lindex $bond 1]]
        set b [lindex $atomtypes $idx]
        if { [string compare $a $b] > 0 } { set t $a; set a $b; set b $t }
        set type [join [list $a $b] "-"]
        lappend newbonds [list [lindex $bond 0] [lindex $bond 1] $type]
    }
    setbondlist $sel type $newbonds
}


# define a new bond or change an existing one.
proc addbond {mol id1 id2 type order} {
    if {$id1 == $id2} {
        vmdcon -err "topo addbond: invalid atom indices: $id1 $id2"
        return
    }

    if {[catch {atomselect $mol "index $id1 $id2"} sel]} {
        vmdcon -err "topo addbond: Invalid atom indices: $sel"
        return
    }

    # make sure we have consistent indexing
    lassign [$sel list] id1 id2

    set bonds [$sel getbonds]
    set bords [$sel getbondorders]
    set btype [$sel getbondtypes]

    set b1 [lindex $bonds 0]
    set b2 [lindex $bonds 1]
    set bo1 [lindex $bords 0]
    set bo2 [lindex $bords 1]
    set bt1 [lindex $btype 0]
    set bt2 [lindex $btype 1]

    # handle the first atom...
    set pos [lsearch -exact -integer $b1 $id2]
    if { $pos < 0} {
        lappend b1 $id2
        lappend bo1 $order
        lappend bt1 $type
    } else {
        set bo1 [lreplace $bo1 $pos $pos $order]
        set bt1 [lreplace $bt1 $pos $pos $type]
    }

    # ...and the second one.
    set pos [lsearch -exact -integer $b2 $id1]
    if { $pos < 0} {
        lappend b2 $id1
        lappend bo2 $order
        lappend bt2 $type
    } else {
        set bo2 [lreplace $bo2 $pos $pos $order]
        set bt2 [lreplace $bt2 $pos $pos $type]
    }

    # and write the modified data back.
    $sel setbonds [list $b1 $b2]
    if {![string equal $order 1.0]} {
        $sel setbondorders [list $bo1 $bo2]
    }
    if {![string equal $type unknown]} {
        $sel setbondtypes [list $bt1 $bt2]
    }
    $sel delete
}

# delete a bond.
proc delbond {mol id1 id2 {type unknown} {order 1.0}} {
    if {[catch {atomselect $mol "index $id1 $id2"} sel]} {
        vmdcon -err "topology delbond: Invalid atom indices: $sel"
        return
    }

    # make sure we have consistent indexing
    lassign [$sel list] id1 id2

    set bonds [$sel getbonds]

    set b1 [lindex $bonds 0]
    set b2 [lindex $bonds 1]

    # handle the first atom...
    set pos [lsearch -exact -integer $b1 $id2]
    if { $pos < 0} {
        ; # bond is not completely within selection. ignore
    } else {
        set b1 [lreplace $b1 $pos $pos]
    }

    # ...and the second one.
    set pos [lsearch -exact -integer $b2 $id1]
    if { $pos < 0} {
        ; # bond is not completely within selection. ignore...
    } else {
        set b2 [lreplace $b2 $pos $pos]
    }

    # and write the modified data back.
    $sel setbonds [list $b1 $b2]
    $sel delete
}

#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009,2010,2011 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topoangles.tcl,v 1.12 2014/08/19 16:45:04 johns Exp $

# return info about angles
# we list and count only angles that are entirely within the selection.
proc angleinfo {infotype sel {flag none}} {

    set numangles 0
    array set angletypes {}
    set atomindex [$sel list]
    set anglelist {}

    foreach angle [join [molinfo [$sel molid] get angles]] {
        lassign $angle t a b c

        if {([lsearch -sorted -integer $atomindex $a] >= 0)          \
                && ([lsearch -sorted -integer $atomindex $b] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $c] >= 0) } {
            set angletypes($t) 1
            incr numangles
            lappend anglelist $angle
        }
    }
    switch $infotype {

        numangles      { return $numangles }
        numangletypes  { return [array size angletypes] }
        angletypenames { return [lsort -ascii [array names angletypes]] }
        getanglelist   { return $anglelist }
        default        { return "bug! shoot the programmer?"}
    }
}

# delete all fully contained angles of the selection.
proc clearangles {sel} {
    set mol [$sel molid]
    set atomindex [$sel list]
    set anglelist {}

    foreach angle [join [molinfo $mol get angles]] {
        lassign $angle t a b c

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0) } {
            lappend anglelist $angle
        }
    }
    molinfo $mol set angles [list $anglelist]
}

# reset angles to data in anglelist
proc setanglelist {sel anglelist} {

    set mol [$sel molid]
    set atomindex [$sel list]
    set newanglelist {}

    # set defaults
    set t unknown; set a -1; set b -1; set c -1

    # preserve all angles definitions that are not contained in $sel
    foreach angle [join [molinfo $mol get angles]] {
        lassign $angle t a b c

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0) } {
            lappend newanglelist $angle
        }
    }

    # append new ones, but only those fully contained in $sel
    foreach angle $anglelist {
        lassign $angle t a b c

        if {([lsearch -sorted -integer $atomindex $a] >= 0)          \
                && ([lsearch -sorted -integer $atomindex $b] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $c] >= 0) } {
            lappend newanglelist $angle
        }
    }

    molinfo $mol set angles [list $newanglelist]
}

# reset angles to data in anglelist
proc retypeangles {sel} {

    set mol [$sel molid]
    set anglelist [angleinfo getanglelist $sel]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newanglelist {}

    foreach angle $anglelist {
        lassign $angle type i1 i2 i3

        set idx [lsearch -sorted -integer $atomindex $i1]
        set a [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i2]
        set b [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i3]
        set c [lindex $atomtypes $idx]

        if { [string compare $a $c] > 0 } { set t $a; set a $c; set c $t }
        set type [join [list $a $b $c] "-"]

        lappend newanglelist [list $type $i1 $i2 $i3]
    }
    setanglelist $sel $newanglelist
}

# reset angles to definitions derived from bonds.
# this includes retyping of the angles.
proc guessangles {sel} {

    set mol [$sel molid]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newanglelist {}

    set bonddata [$sel getbonds]

    # preserve all angles definitions that are not fully contained in $sel
    foreach angle [angleinfo getanglelist $sel] {
        lassign $angle t a b c

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0) } {
            lappend newanglelist $angle
        }
    }

    # a topological angle is defined by two bonds that share an atom
    # bound to it that are not the bond itself
    foreach bonds $bonddata aidx $atomindex atyp $atomtypes {
        set nbnd [llength $bonds]
        for {set i 0} {$i < $nbnd-1} {incr i} {
            for {set j [expr {$i+1}]} {$j < $nbnd} {incr j} {
                set b1idx [lindex $bonds $i]
                set idx [lsearch -sorted -integer $atomindex $b1idx]
                set b1typ [lindex $atomtypes $idx]
                set b2idx [lindex $bonds $j]
                set idx [lsearch -sorted -integer $atomindex $b2idx]
                set b2typ [lindex $atomtypes $idx]
                if { ([string compare $b1typ $b2typ] > 0) } {
                    set t1 $b1typ; set b1typ $b2typ; set b2typ $t1
                    set t2 $b1idx; set b1idx $b2idx; set b2idx $t2
                }
                set type [join [list $b1typ $atyp $b2typ] "-"]

                # append only angles that are full contained in $sel
                if {([lsearch -sorted -integer $atomindex $b1idx] >= 0)          \
                        && ([lsearch -sorted -integer $atomindex $aidx] >= 0)   \
                        && ([lsearch -sorted -integer $atomindex $b2idx] >= 0) } {
                    lappend newanglelist [list $type $b1idx $aidx $b2idx]
                }
            }
        }
    }
    molinfo $mol set angles [list $newanglelist]
}

# define a new angle or change an existing one.
proc addangle {mol id1 id2 id3 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3"} sel]} {
        vmdcon -err "topology addangle: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id1 > $id3} {set t $id1 ; set id1 $id3 ; set id3 $t }

    set angles [join [molinfo $mol get angles]]
    lappend angles [list $type $id1 $id2 $id3]
    $sel delete
    molinfo $mol set angles [list $angles]
}

# delete an angle.
proc delangle {mol id1 id2 id3 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3"} sel]} {
        vmdcon -err "topology delangle: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id1 > $id3} {set t $id1 ; set id1 $id3 ; set id3 $t }

    set newanglelist {}
    foreach angle [join [molinfo $mol get angles]] {
        lassign $angle t a b c
        if { ($a != $id1) || ($b != $id2) || ($c != $id3) } {
            lappend newanglelist $angle
        }
    }
    $sel delete
    molinfo $mol set angles [list $newanglelist]
}

#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009,2010,2011 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topodihedrals.tcl,v 1.11 2014/08/19 16:45:04 johns Exp $

# return info about dihedrals
# we list and count only dihedrals that are entirely within the selection.
proc dihedralinfo {infotype sel {flag none}} {

    set numdihedrals 0
    array set dihedraltypes {}
    set atomindex [$sel list]
    set dihedrallist {}

    foreach dihedral [join [molinfo [$sel molid] get dihedrals]] {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atomindex $a] >= 0)          \
                && ([lsearch -sorted -integer $atomindex $b] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $c] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $d] >= 0) } {
            set dihedraltypes($t) 1
            incr numdihedrals
            lappend dihedrallist $dihedral
        }
    }
    switch $infotype {

        numdihedrals      { return $numdihedrals }
        numdihedraltypes  { return [array size dihedraltypes] }
        dihedraltypenames { return [lsort -ascii [array names dihedraltypes]] }
        getdihedrallist   { return $dihedrallist }
        default        { return "bug! shoot the programmer?"}
    }
}

# delete all contained dihedrals of the selection.
proc cleardihedrals {sel} {
    set mol [$sel molid]
    set atomindex [$sel list]
    set dihedrallist {}

    foreach dihedral [join [molinfo $mol get dihedrals]] {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0)   \
                || ([lsearch -sorted -integer $atomindex $d] < 0) } {
            lappend dihedrallist $dihedral
        }
    }
    molinfo $mol set dihedrals [list $dihedrallist]
}

# reset dihedrals to data in dihedrallist
proc setdihedrallist {sel dihedrallist} {

    set mol [$sel molid]
    set atomindex [$sel list]
    set newdihedrallist {}

    # set defaults
    set t unknown; set a -1; set b -1; set c -1; set d -1

    # preserve all dihedrals definitions that are not fully contained in $sel
    foreach dihedral [join [molinfo $mol get dihedrals]] {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0)   \
                || ([lsearch -sorted -integer $atomindex $d] < 0) } {
            lappend newdihedrallist $dihedral
        }
    }

    # append new ones, but only those contained in $sel
    foreach dihedral $dihedrallist {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atomindex $a] >= 0)          \
                && ([lsearch -sorted -integer $atomindex $b] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $c] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $d] >= 0) } {
            lappend newdihedrallist $dihedral
        }
    }

    molinfo $mol set dihedrals [list $newdihedrallist]
}

# reset dihedrals to data in dihedrallist
proc retypedihedrals {sel} {

    set mol [$sel molid]
    set dihedrallist [dihedralinfo getdihedrallist $sel]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newdihedrallist {}

    foreach dihedral $dihedrallist {
        lassign $dihedral type i1 i2 i3 i4

        set idx [lsearch -sorted -integer $atomindex $i1]
        set a [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i2]
        set b [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i3]
        set c [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i4]
        set d [lindex $atomtypes $idx]

        if { ([string compare $b $c] > 0) \
                 || ( [string equal $b $c] && [string compare $a $d] > 0 ) } {
            set t $a; set a $d; set d $t
            set t $b; set b $c; set c $t
            set t $i1; set i1 $i4; set i4 $t
            set t $i2; set i2 $i3; set i3 $t
        }
        set type [join [list $a $b $c $d] "-"]

        lappend newdihedrallist [list $type $i1 $i2 $i3 $i4]
    }
    setdihedrallist $sel $newdihedrallist
}


# reset dihedrals to definitions derived from bonds.
# this includes retyping of the dihedrals.
proc guessdihedrals {sel} {

    set mol [$sel molid]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newdihedrallist {}

    set bondlist [bondinfo getbondlist $sel]
    set bonddata [$sel getbonds]

    # preserve all dihedrals definitions that are not fully contained in $sel
    foreach dihedral [join [molinfo $mol get dihedrals]] {
        lassign $dihedral t a b c d

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0)   \
                || ([lsearch -sorted -integer $atomindex $d] < 0) } {
            lappend newdihedrallist $dihedral
        }
    }

    # a topological dihedral is defined by a bond and atoms
    # bound to it that are not the bond itself
    foreach bond $bondlist {
        lassign $bond b1 b2
        set b1idx [lsearch -sorted -integer $atomindex $b1]
        set b1typ [lindex $atomtypes $b1idx]
        set b2idx [lsearch -sorted -integer $atomindex $b2]
        set b2typ [lindex $atomtypes $b2idx]
        foreach o1 [lindex $bonddata $b1idx] {
            foreach o2 [lindex $bonddata $b2idx] {
                if {($o1 == $b1) || ($o2 == $b1) || ($o1 == $b2) || ($o2 == $b2)} {
                    continue
                }
                set o1idx [lsearch -sorted -integer $atomindex $o1]
                set o1typ [lindex $atomtypes $o1idx]
                set o2idx [lsearch -sorted -integer $atomindex $o2]
                set o2typ [lindex $atomtypes $o2idx]
                if { ([string compare $b1typ $b2typ] > 0) \
                 || ( [string equal $b1typ $b2typ]
                      && [string compare $o1typ $o2typ] > 0 ) } {
                    set type [join [list $o2typ $b2typ $b1typ $o1typ] "-"]
                    lappend newdihedrallist [list $type $o2 $b2 $b1 $o1]
                } else {
                    set type [join [list $o1typ $b1typ $b2typ $o2typ] "-"]
                    lappend newdihedrallist [list $type $o1 $b1 $b2 $o2]
                }
            }
        }
    }
    setdihedrallist $sel $newdihedrallist
}


# define a new dihedral or change an existing one.
proc adddihedral {mol id1 id2 id3 id4 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3 $id4"} sel]} {
        vmdcon -err "topology adddihedral: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id2 > $id3} {
        set t $id2 ; set id2 $id3 ; set id3 $t
        set t $id1 ; set id1 $id4 ; set id4 $t
    }

    set dihedrals [join [molinfo $mol get dihedrals]]
    lappend dihedrals [list $type $id1 $id2 $id3 $id4]
    molinfo $mol set dihedrals [list $dihedrals]
    # this is not (yet) required
    $sel delete
    return
}

# delete a dihedral.
proc deldihedral {mol id1 id2 id3 id4 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3 $id4"} sel]} {
        vmdcon -err "topology deldihedral: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id2 > $id3} {
        set t $id2 ; set id2 $id3 ; set id3 $t
        set t $id1 ; set id1 $id4 ; set id4 $t
    }

    set newdihedrallist {}
    foreach dihedral [join [molinfo $mol get dihedrals]] {
        lassign $dihedral t a b c d
        if { ($a != $id1) || ($b != $id2) || ($c != $id3) || ($d != $id4) } {
            lappend newdihedrallist $dihedral
        }
    }
    molinfo $mol set dihedrals [list $newdihedrallist]
    # this is not (yet) required
    $sel delete
    return
}

#!/usr/bin/tclsh
# This file is part of TopoTools, a VMD package to simplify
# manipulating bonds other topology related properties.
#
# Copyright (c) 2009,2010,2011 by Axel Kohlmeyer <akohlmey@gmail.com>
# $Id: topoimpropers.tcl,v 1.11 2014/08/19 16:45:04 johns Exp $

# return info about impropers
# we list and count only impropers that are entirely within the selection.
proc improperinfo {infotype sel {flag none}} {

    set numimpropers 0
    array set impropertypes {}
    set atomindex [$sel list]
    set improperlist {}

    foreach improper [join [molinfo [$sel molid] get impropers]] {
        lassign $improper t a b c d

        if {([lsearch -sorted -integer $atomindex $a] >= 0)          \
                && ([lsearch -sorted -integer $atomindex $b] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $c] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $d] >= 0) } {
            set impropertypes($t) 1
            incr numimpropers
            lappend improperlist $improper
        }
    }
    switch $infotype {

        numimpropers      { return $numimpropers }
        numimpropertypes  { return [array size impropertypes] }
        impropertypenames { return [lsort -ascii [array names impropertypes]] }
        getimproperlist   { return $improperlist }
        default        { return "bug! shoot the programmer?"}
    }
}

# delete all contained impropers of the selection.
proc clearimpropers {sel} {
    set mol [$sel molid]
    set atomindex [$sel list]
    set improperlist {}

    foreach improper [join [molinfo $mol get impropers]] {
        lassign $improper t a b c d

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0)   \
                || ([lsearch -sorted -integer $atomindex $d] < 0) } {
            lappend improperlist $improper
        }
    }
    molinfo $mol set impropers [list $improperlist]
}

# reset impropers to data in improperlist
proc setimproperlist {sel improperlist} {

    set mol [$sel molid]
    set atomindex [$sel list]
    set newimproperlist {}

    # set defaults
    set t unknown; set a -1; set b -1; set c -1; set d -1

    # preserve all impropers definitions that are not contained in $sel
    foreach improper [join [molinfo $mol get impropers]] {
        lassign $improper t a b c d

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0)   \
                || ([lsearch -sorted -integer $atomindex $d] < 0) } {
            lappend newimproperlist $improper
        }
    }

    # append new ones, but only those contained in $sel
    foreach improper $improperlist {
        lassign $improper t a b c d

        if {([lsearch -sorted -integer $atomindex $a] >= 0)          \
                && ([lsearch -sorted -integer $atomindex $b] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $c] >= 0)   \
                && ([lsearch -sorted -integer $atomindex $d] >= 0) } {
            lappend newimproperlist $improper
        }
    }

    molinfo $mol set impropers [list $newimproperlist]
}

# reset impropers to data in improperlist
proc retypeimpropers {sel} {

    set mol [$sel molid]
    set improperlist [improperinfo getimproperlist $sel]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newimproperlist {}

    foreach improper $improperlist {
        lassign $improper type i1 i2 i3 i4

        set idx [lsearch -sorted -integer $atomindex $i1]
        set a [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i2]
        set b [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i3]
        set c [lindex $atomtypes $idx]
        set idx [lsearch -sorted -integer $atomindex $i4]
        set d [lindex $atomtypes $idx]

        if { ([string compare $b $c] > 0) \
                 || ( [string equal $b $c] && [string compare $a $d] > 0 ) } {
            set t $a; set a $d; set d $t
            set t $b; set b $c; set c $t
            set t $i1; set i1 $i4; set i4 $t
            set t $i2; set i2 $i3; set i3 $t
        }
        set type [join [list $a $b $c $d] "-"]

        lappend newimproperlist [list $type $i1 $i2 $i3 $i4]
    }
    setimproperlist $sel $newimproperlist
}

# reset impropers to definitions derived from bonds.
# this includes retyping of the impropers.
# this step is different from guessing angles or dihedrals,
# as we are only looking for definitions that are unusual.

proc guessimpropers {sel {flags {}}} {
    # default tolerance is 5 degrees from planar
    set tolerance 5

    # parse optional flags
    foreach {key value} $flags {
        switch -- $key {
            tol -
            tolerance {set tolerance  $value}
            default {
                vmdcon -err "guessimpropers: unknown flag: $key"
                return -1
            }
        }
    }

    set mol [$sel molid]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newimproperlist {}

    set bonddata [$sel getbonds]
    set minangle [expr {180.0 - $tolerance}]

    # preserve all impropers definitions that are not fully contained in $sel
    foreach improper [join [molinfo $mol get impropers]] {
        lassign $improper t a b c d

        if {([lsearch -sorted -integer $atomindex $a] < 0)          \
                || ([lsearch -sorted -integer $atomindex $b] < 0)   \
                || ([lsearch -sorted -integer $atomindex $c] < 0)   \
                || ([lsearch -sorted -integer $atomindex $d] < 0) } {
            lappend newimproperlist $improper
        }
    }

    # a topological improper is defined by three bonds connected to
    # the same atom and their dihedral being almost in plane.
    foreach bonds $bonddata aidx $atomindex atyp $atomtypes {
        set nbnd [llength $bonds]
        if {$nbnd == 3} {
            lassign $bonds b1 b2 b3
            set ang [expr {abs([measure imprp [list $b1 $b2 $aidx $b3] molid $mol])}]
            if {$ang > $minangle} {
                set b1idx [lsearch -sorted -integer $atomindex $b1]
                set b1typ [lindex $atomtypes $b1idx]
                set b2idx [lsearch -sorted -integer $atomindex $b2]
                set b2typ [lindex $atomtypes $b2idx]
                set b3idx [lsearch -sorted -integer $atomindex $b3]
                set b3typ [lindex $atomtypes $b3idx]

                if {([string compare $b1typ $b2typ]) > 0} {
                    set t1 $b1typ; set b1typ $b2typ; set b2typ $t1
                    set t2 $b1; set b1 $b2; set b2 $t2
                }
                if {([string compare $b2typ $b3typ]) > 0} {
                    set t1 $b2typ; set b2typ $b3typ; set b3typ $t1
                    set t2 $b2; set b2 $b3; set b3 $t2
                }
                if {([string compare $b1typ $b2typ]) > 0} {
                    set t1 $b1typ; set b1typ $b2typ; set b2typ $t1
                    set t2 $b1; set b1 $b2; set b2 $t2
                }
                set type [join [list $b1typ $b2typ $atyp $b3typ] "-"]
                lappend newimproperlist [list $type $b1 $b2 $aidx $b3]
            }
        }
    }
    setimproperlist $sel $newimproperlist
}

# define a new improper or change an existing one.
proc addimproper {mol id1 id2 id3 id4 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3 $id4"} sel]} {
        vmdcon -err "topology addimproper: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id2 > $id3} {
        set t $id2 ; set id2 $id3 ; set id3 $t
        set t $id1 ; set id1 $id4 ; set id4 $t
    }

    set impropers [join [molinfo $mol get impropers]]
    lappend impropers [list $type $id1 $id2 $id3 $id4]
    $sel delete
    molinfo $mol set impropers [list $impropers]
}

# delete a improper.
proc delimproper {mol id1 id2 id3 id4 {type unknown}} {
    if {[catch {atomselect $mol "index $id1 $id2 $id3 $id4"} sel]} {
        vmdcon -err "topology delimproper: Invalid atom indices: $sel"
        return
    }

    # canonicalize indices
    if {$id2 > $id3} {
        set t $id2 ; set id2 $id3 ; set id3 $t
        set t $id1 ; set id1 $id4 ; set id4 $t
    }

    set newimproperlist {}
    foreach improper [join [molinfo $mol get impropers]] {
        lassign $improper t a b c d
        if { ($a != $id1) || ($b != $id2) || ($c != $id3) || ($d != $id4) } {
            lappend newimproperlist $improper
        }
    }
    $sel delete
    molinfo $mol set impropers [list $newimproperlist]
}


set datafile [format "%s.data" $::argv]
set psffile [format "%s.psf" $::argv]
readlammpsdata $datafile full
animate write psf $psffile
exit

