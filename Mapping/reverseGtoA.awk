BEGIN { 
    FS="\t"
    OFS="\t"
} 

/^@/ {print}

$0!~/^@/ {
    if ( $1 == old_name ) {
        print old_string
        next
    }
    
    while ( line2[1] != $1 ) {
        getline line < fbed
        split(line, line2, "\t")
    }

    if ( line2[2] != "NA" ) {
        n_sites = split(line2[2], sites, "_")
        n_base = split($10, read_a, "")
        for (i=1;i<=n_sites;i++) {
            read_a[sites[i]+1] = "A"
        }
    
        read_s = read_a[1]
        for (i=2; i<=n_base; i++) 
            read_s = read_s read_a[i]
        $10 = read_s
    }
    
    print $0
    old_name = $1
    old_string = $0
}
