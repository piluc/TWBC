include("TemporalWalkBetweenness.jl")

function read_tsb_values(fn)
    f = open("local/sea/scores/" * fn * "/tsb.txt", "r")
    l = readlines(f)
    tsb = zeros(length(l))
    for i in 1:lastindex(l)
        tsb[i] = parse(Float64, l[i])
    end
    return tsb
end

function compare_tsb_twbc(fn)
    println(fn)
    tsb = read_tsb_values(fn)
    start_time = time()
    twbc = temporal_walk_betweenness("local/patg/" * fn * ".patg", " ", false)
    exec_time = round(time() - start_time; digits=4)
    logging("TWBC computed  in " * string(exec_time) * " seconds", true, false)
    for i in 1:lastindex(tsb)
        if (!isapprox(tsb[i], twbc[i]; atol=0.00001))
            logging(string(tsb[i]) * " != " * string(twbc[i]) * " at " * string(i), true, false)
            return
        end
    end
    logging("TSB = TWBC", true, false)
end

function compare_tsb_twbc()
    # nn = ["001_highschool_2011", "002_highschool_2012", "003_highschool_2013", "004_hospital_ward", "005_hypertext_2009", "006_primary_school", "01_venice", "02_college_msg", "03_email_eu", "04_bordeaux", "05_infectious", "06_topology", "07_wiki_elections", "08_digg_reply", "09_slashdot_reply", "10_facebook_wall", "12_SMS"]
    nn = ["05_infectious", "06_topology", "07_wiki_elections", "08_digg_reply", "09_slashdot_reply", "10_facebook_wall", "12_SMS"]
    for ni in 1:lastindex(nn)
        compare_tsb_twbc(nn[ni])
    end
end
