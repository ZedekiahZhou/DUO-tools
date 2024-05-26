BEGIN { FS=" " } 

NR % 4 == 1 {  # only keep seq ID, delete "@"
    print $1 >> changefastq
    sub("@", "", $1)
    printf "%s\t", $1
}

NR % 4 == 2 {
    len=split($1,arr,"")
    j=0
    for ( i=1;i<=len;i++ ) {  # record A sites to arr
        if (arr[i]=="A") {
            j+=1
            pos[j]=i
        }
    }

    if ( j == 0) 
        print "NA"
    else {
        seq=pos[1]-1
        for (i=2;i<=j;i++) {  # concat the positions with '_'
            seq=seq "_" pos[i]-1
        }
        print seq
    }
    gsub("A", "G", $1)
    print $1 >> changefastq
}

NR % 4 == 3 || NR % 4 == 0 { print $0 >> changefastq }