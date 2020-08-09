using LasIO
using LazIO
using FileIO
using StaticArrays
using Plots
using JSON

## ------------FUNCTIONS-------------- ##
"Group points of pulse, dictinary index --> tag,intensity,z value, return number"
function grouping(points,header)
    prev_point = points[1]
    group_pulse = Dict{Int,Tuple{String,UInt16,Float64,UInt8}}() #dictionary with indexes and tags
    n = 0
    for (i,p) in enumerate(points)
        inten = intensity(p)
        z = zcoord(p,header)
        rn = return_number(p)
        nr = number_of_returns(p)
        dt = gps_time(p) - gps_time(prev_point)
        if (rn == n ) && (nr > 1) #for every pulse
            if rn == 1
                group_pulse[i] = ("firstpoint",inten,z,rn)
            elseif rn < nr && dt < 1e-7
                group_pulse[i] = ("nextpoint",inten,z,rn)
            elseif rn == nr && dt < 1e-7
                group_pulse[i] = ("lastpoint",inten,z,rn)
            end
            if nr == rn
                n = 0
            end
        elseif nr == 1 && rn == 1 #cases that I have JUST one point!
            group_pulse[i] = ("firstpoint",inten,z,rn)
            n = 0
        elseif nr > 1 && rn == 1 #cases that I have JUST one point!
            group_pulse[i] = ("firstpoint",inten,z,rn)
            n = 0
        elseif nr == rn && nr!=1 #cases that I have remaining points like rn=3, nr=3 from other pulse!
            group_pulse[i] = ("lastpoint",inten,z,rn)
            n = 0
        else
            group_pulse[i] = ("unclassified",inten,z,rn)
            n = 0
        end
        n+=1
        prev_point = p
    end
    group_pulse
end

"Write LAS file"
function write_output(writefiles, ds, points)
    LazIO.write(writefiles, ds.header) do io
        for (i,p) in enumerate(ds)
            if i in points
                p.classification = UInt(5)
                LazIO.writepoint(io,p)
            end
        end
    end
end


## ---------------MAIN-----------------##
"read - import dataset from a json file"
workdir = dirname(@__FILE__)
open(joinpath(workdir,"pulse.json"), "r") do f
    global mydict = Dict(JSON.parse(read(f,String)))  # parse and transform data
end

#get the filename
filenn_NL1 = joinpath(workdir,mydict["dataname"] * ".laz")
filenn_NL1_las = File{format"LAZ_"}(filenn_NL1)

#read specific LAS file
header, points = LazIO.load(filenn_NL1_las) #point cloud as LAS points format
n = length(points) #number of points

#read dataset
ds = LazIO.open(filenn_NL1) #dataset

#export plots
# path = "C:/Users/alexandr/OneDrive - Stichting Deltares/Desktop/Thesis/P3_processing/with_greenstowa/Working_to/Plots_378_NL3"
# println("Export path: Ok")
# name = "141_NL3"

#Grouping_points per pulse into 3 and 2 points
group_pulses = grouping(points,header) #point id --> tag, intensity, z_val, return number (rn)
allkeys = collect(keys(group_pulses)) #keys
p3_pulse = Vector{Vector{Float64}}[] ; p2_pulse = Vector{Vector{Float64}}[] #store triplets  and doubles of pulses
bottom_points_3 = Int64[]; bottom_points_2 = Int64[]
all_bottoms = Int64[]
for (i,key) in enumerate(allkeys)
    lex = group_pulses[key]
    if lex[1] == "nextpoint"
        if group_pulses[allkeys[i-1]][1] == "firstpoint" && group_pulses[allkeys[i+1]][1] == "lastpoint" && group_pulses[allkeys[i+1]][4] != group_pulses[allkeys[i]][4] #check the return number of last and second point

            #intensities
            sec_int = group_pulses[allkeys[i]][2]
            fir_int = group_pulses[allkeys[i-1]][2]
            last_int = group_pulses[allkeys[i+1]][2]

            #Check if LAST point's intensity >= MIDDLE point's intensity
            if (last_int/sec_int) >= 1
                push!(bottom_points_3,key)
                push!(all_bottoms,key)
            end

            #z_values
            sec_z = group_pulses[allkeys[i]][3]
            fir_z = group_pulses[allkeys[i-1]][3]
            last_z = group_pulses[allkeys[i+1]][3]

            #return number (rn)
            sec_rn = group_pulses[allkeys[i]][4]
            fir_rn = group_pulses[allkeys[i-1]][4]
            last_rn = group_pulses[allkeys[i+1]][4]

            push!(p3_pulse,[[fir_int,fir_z,fir_rn], [sec_int,sec_z,sec_rn], [last_int,last_z,last_rn]])
        end
    elseif lex[1] == "lastpoint"
        if group_pulses[allkeys[i-1]][1] == "firstpoint"

            #intensities
            fir_int = group_pulses[allkeys[i-1]][2]
            sec_int = group_pulses[allkeys[i]][2]

            #z_values
            fir_z = group_pulses[allkeys[i-1]][3]
            sec_z = group_pulses[allkeys[i]][3]

            #return number (rn)
            fir_rn = group_pulses[allkeys[i-1]][4]
            sec_rn = group_pulses[allkeys[i]][4]

            push!(p2_pulse,[[fir_int,fir_z,fir_rn], [sec_int,sec_z,sec_rn]])
            push!(bottom_points_2,key)
            push!(all_bottoms,key)
        end
    end
end

println("length of p3 pulses: ",length(p3_pulse))
println("length of p2 pulses: ",length(p2_pulse))

#write output files
writefiles = joinpath(workdir, mydict["dataname"] * mydict["outputname"] * ".laz")
write_output(writefiles, ds, all_bottoms)
