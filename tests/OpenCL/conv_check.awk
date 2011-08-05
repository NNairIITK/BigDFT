BEGIN {
    max_discrepancy = 2.e-10
    print max_discrepancy
    discrepancy = 0.0
    end = "\033[m"
    maxline = 155
    }
/nbelem/ {
    if ($6 > discrepancy) {
        discrepancy = $6
    }
    if (discrepancy > max_discrepancy) {
        start = "\033[0;31m"
        print $0;
        printf("%sMax discrepancy : %7.1e > (%7.1e)    failed%s\n",start,discrepancy,max_discrepancy,end)
        exit 1
    }
}

END {
        if (NR != maxline) {
            start = "\033[0;31m"
            printf("The number of line is not correct (%d /= %d)\n",NR,maxline);
            printf("%sMax discrepancy : nan > (%7.1e)    failed%s\n",start,max_discrepancy,end)
            exit 1
            }
        start = "\033[0;32m"
        printf("%sMax discrepancy : %7.1e < (%7.1e) succeeded%s\n",start,discrepancy,max_discrepancy,end)
}
