package require pbctools

proc read_xyz { fname } {
    mol new $fname waitfor all
    set molid [molinfo top get id]
    set fp [open $fname "r"]
    set data [read $fp]
    close $fp
    set lines [split $data \n]
    set boxdat {}
    foreach line $lines {
        set dat {}
        foreach l $line { lappend dat $l }
        if { [string equal [lindex $dat 0] "BOX"]} {
            lappend boxdat [lrange $dat 1 end]
        }
    }
    for { set f 0 } { $f < [molinfo top get numframes] } { incr i } {
        pbc set [lindex $boxdat $f] -first $f -last $f
    }
    return $molid
}
