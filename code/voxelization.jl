using LasIO
using LazIO
using FileIO
using StaticArrays
using GeoArrays; const G = GeoArrays
using Colors
using PyPlot
using PyCall
using Plots
using ArchGDAL ; const AG = ArchGDAL ## make AG the shortcut to ArchGDAL, by doing this I call AG. rather than ArchGDAL.
using JSON

## ------------FUNCTIONS-------------- ##
"Store coordinates"
function coordinates(points,header)
    # coords_array = SArray{Tuple{4},Float64,1,4}[]
    coords_dict = Dict{Int64,Array{Float64,1}}()
    for (ind,p) in enumerate(points)
        x = xcoord(p,header)
        y = ycoord(p,header)
        z = zcoord(p,header)
        inten = intensity(p)
        ret = return_number(p)
        nret = number_of_returns(p)
        classificc = classification(p)
        gps = gps_time(p)
        coords_dict[ind] = [x,y,z,inten]
        # push!(coords_array,(ind, x, y, z))
    end
    return coords_dict
end

"Returns minimum boundary point"
function minboundrp(point::SArray{Tuple{3},Float64,1,3},voxelsize::SArray{Tuple{3},Float64,1,3})
    x,y,z = point
    zpx = x - voxelsize[1]
    zpy = y - voxelsize[2]
    zpz = z - voxelsize[3]
    newpoint = SArray{Tuple{3},Float64,1,3}(zpx,zpy,zpz)
    return newpoint
end

"Returns maximum boundary point"
function maxboundrp(point::SArray{Tuple{3},Float64,1,3},voxelsize::SArray{Tuple{3},Float64,1,3})
    x,y,z = point
    zpx = x + voxelsize[1]
    zpy = y + voxelsize[2]
    zpz = z + voxelsize[3]
    newpoint = SArray{Tuple{3},Float64,1,3}(zpx,zpy,zpz)
    return newpoint
end

"Local maxima, returns indexes"
function findlocalmaxima_indexes(signal)
    inds = Int[]
    if length(signal)>1
        if signal[1]>signal[2]
            push!(inds,1)
        end
        for i=2:length(signal)-1
            if signal[i-1]<signal[i]>signal[i+1]
                push!(inds,i)
            end
        end
        if signal[end]>signal[end-1]
            push!(inds,length(signal))
        end
    end
    inds
end

"Creates a histogram for every voxel"
function myhistogram(gr_points_to_voxels)
    alldepth = Float64[]; count_voxs = Int64[]
    allkeys = collect(keys(gr_points_to_voxels))
    for i in allkeys
        voxel = gr_points_to_voxels[i]
        len = length(voxel)
        # depth = Float64[] #vector with Floats
        for i in voxel
            z = coords_dict[i][3]
            push!(alldepth,z)
        end
        push!(count_voxs,len)
    end
    len = length(alldepth)
    minn = minimum(alldepth)
    maxx = maximum(alldepth)
    return minn,maxx,len,alldepth, allkeys
end

"Range of the z_start and z_stop of all the histogram bins of a voxel"
function ranges(rang1,rang2,voxel)
    ids =  Dict{Int64,Float64}()
    for i in voxel
        if rang1 <= coords_dict[i][3] <= rang2
            ids[i] = coords_dict[i][3]
        end
    end
    return ids
end

"Write confindence values as an attribute of a point"
function write_confidence(writeoutput, ds, confiden_vals, median_nden, median_nint, mean_ndis)
    confidenkeys = collect(keys(confiden_vals))
    LazIO.write(writeoutput, ds.header) do io
        for (i,p) in enumerate(ds) #I need to interate only to dataset, since it is Lazpoints
            if i in confidenkeys
                #1. density, 2. distance, 3. intensity
                if confiden_vals[i][1] > median_nden &&   confiden_vals[i][2] > mean_ndis && confiden_vals[i][3] > median_nint
                    p.classification = UInt(1)
                elseif confiden_vals[i][1] > median_nden && confiden_vals[i][2] <= mean_ndis  && confiden_vals[i][3] > median_nint
                    p.classification = UInt(2)
                elseif confiden_vals[i][1] > median_nden &&   confiden_vals[i][2] > mean_ndis && confiden_vals[i][3] <= median_nint
                    p.classification = UInt(3)
                elseif confiden_vals[i][1] <= median_nden &&  confiden_vals[i][2] > mean_ndis &&  confiden_vals[i][3] > median_nint
                    p.classification = UInt(4)
                elseif  confiden_vals[i][1] <= median_nden  &&  confiden_vals[i][2] > mean_ndis &&  confiden_vals[i][3] <= median_nint
                    p.classification = UInt(5)
                elseif confiden_vals[i][1] > median_nden && confiden_vals[i][2] <= mean_ndis  && confiden_vals[i][3] <= median_nint
                    p.classification = UInt(6)
                elseif confiden_vals[i][1] <= median_nden &&  confiden_vals[i][2] <= mean_ndis  &&confiden_vals[i][3] > median_nint
                    p.classification = UInt(7)
                elseif  confiden_vals[i][1] <= median_nden  &&  confiden_vals[i][2] <= mean_ndis   && confiden_vals[i][3] <= median_nint
                    p.classification = UInt(8)
                end
                LazIO.writepoint(io,p)
            end
        end
    end
end

function confidence(ds, confiden_vals, median_nden, median_nint, mean_ndis)
    confidenkeys = collect(keys(confiden_vals))
    conf = Dict()
    for (i,p) in enumerate(ds) #I need to interate only to dataset, since it is Lazpoints
        if i in confidenkeys
            #1. density, 2. distance, 3. intensity
            if confiden_vals[i][1] > median_nden &&   confiden_vals[i][2] > mean_ndis && confiden_vals[i][3] > median_nint
                confid = 1
            elseif confiden_vals[i][1] > median_nden && confiden_vals[i][2] <= mean_ndis  && confiden_vals[i][3] > median_nint
                confid = 2
            elseif confiden_vals[i][1] > median_nden &&   confiden_vals[i][2] > mean_ndis && confiden_vals[i][3] <= median_nint
                confid = 3
            elseif confiden_vals[i][1] <= median_nden &&  confiden_vals[i][2] > mean_ndis &&  confiden_vals[i][3] > median_nint
                confid = 4
            elseif  confiden_vals[i][1] <= median_nden  &&  confiden_vals[i][2] > mean_ndis &&  confiden_vals[i][3] <= median_nint
                confid = 5
            elseif confiden_vals[i][1] > median_nden && confiden_vals[i][2] <= mean_ndis  && confiden_vals[i][3] <= median_nint
                confid = 6
            elseif confiden_vals[i][1] <= median_nden &&  confiden_vals[i][2] <= mean_ndis  &&confiden_vals[i][3] > median_nint
                confid = 7
            elseif  confiden_vals[i][1] <= median_nden  &&  confiden_vals[i][2] <= mean_ndis   && confiden_vals[i][3] <= median_nint
                confid = 8
            end
            conf[i] = [confid]
        end
    end
    return conf
end






## ---------------MAIN-----------------##
"read - import dataset from a json file"
workdir = dirname(@__FILE__)
open(joinpath(workdir,"voxelization.json"), "r") do f
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
bb = LasIO.boundingbox(header) #boundary box of the pointcloud
coords_dict = coordinates(points,header) #coordinates


#### "VOXELIZATION_ALGORITHM, Nourian et.al." ####
#voxelsize in x,y dimensions
minb = bb[3]; maxb = bb[6]
voxel_z = round((maxb - minb),digits=1)
voxelsize = SArray{Tuple{3},Float64,1,3}(mydict["vx"], mydict["vy"], voxel_z)

#new boundary - lower left and upper right
bmin = SArray{Tuple{3},Float64,1,3}(bb.xmin, bb.ymin, bb.zmin)
bmax = SArray{Tuple{3},Float64,1,3}(bb.xmax, bb.ymax, bb.zmax)
nb_min = minboundrp(bmin,voxelsize) #new boundary- lower left
nb_max = maxboundrp(bmax,voxelsize) #new boundary- upper right

#dictionary all the points (indexes) --> i,j,k
vox_dict = Dict{Int64,Array{Int64,1}}()
for (ind,p) in enumerate(coords_dict)
    i = round(abs(coords_dict[ind][1] - nb_min[1])/voxelsize[1])
    j = round(abs(coords_dict[ind][2] - nb_min[2])/voxelsize[2])
    k = round(abs(coords_dict[ind][3] - nb_min[3])/voxelsize[3])
    vox_dict[ind] = [i,j,k]
end

#UNIQUE VALUES in the dictionary
#Store INDEXES of the points that correspond to EVERY Voxel
val = unique(values(vox_dict))
#dictionary all the voxels (keys) --> all the indexes of the points PER Voxel (values)
gr_points_to_voxels = Dict{Int64,Array{Int64,1}}()
for (i,every) in enumerate(val)
    result = filter((k,v)->v == every,vox_dict)
    gr_points_to_voxels[i] = collect(keys(result))
end

#PLOTS
#Return some characteristics of a histogram
minn, maxx, len, alldepth, allkeys = myhistogram(gr_points_to_voxels)

nplot = length(gr_points_to_voxels) # number of plots
num_plots = [i for i in range(1,stop = length(gr_points_to_voxels)-(length(gr_points_to_voxels)-nplot))] #selecting specific number of plots for plotting

#Create histogram plots and save the required points
c = 1
imin =  Dict{Int64,Array{Float64,1}}() # store the point id --> density, distance, intensity
watersurface =  Int64[] # just the watersurface points
for i in allkeys
    if c <= length(num_plots)
        voxel = gr_points_to_voxels[i] #points of a voxel
        npoints = length(voxel) #amount of points per voxel
        if length(voxel) > 1 #if there are points
            zvals = Float64[] #vector with Floats
            for i in voxel
                push!(zvals,coords_dict[i][3])
            end
            len = length(zvals); minn = minimum(zvals); maxx = maximum(zvals)
            dis = abs(maxx - minn);  width = 0.05  #5cm spacing
            nb = Int64(round(dis/width)) #number of bins per voxel

            if nb > 1
                #Create a histogram PER Voxel, with a specific Number of Bins (nb)
                plt.figure()
                plt.hist(zvals,bins=nb) #check it here
                plt.ylabel("Count of points")
                plt.xlabel("Z value (m)")
                plt.title("Voxel "*string(i)*": Count of points with respect to z value")

                #Extract info from a histogram
                counts, bins, bars = plt.hist(zvals, bins=nb)
                bins = [round(i,digits=2) for i in bins] #round them to 2 digits
                println("mycounts : ",counts)
                println("my bins : ",bins)
                println("Voxel: ",i)
                println()

                #1. MAX BIN
                # find all INDEXES the peaks from counts
                peaks = findlocalmaxima_indexes(counts)
                println("My peaks in the histogram: ",peaks)
                # just one peak
                if length(peaks) == 1
                    val = peaks[1] #extract that element

                    density = counts[val]/npoints #density
                    distance = -1.0 #distance

                    rang = bins[val:val+1] #find the values of Z range of that bin

                    if bins[val] >= mean(bins) #if the HIGHEST peak is on the RIGHT side of the Mean

                        # find the lowest point of the BIN
                        ids_above = ranges(rang[1],rang[2],voxel)
                        if length(ids_above) > 0
                            minpoint_high = findmin(ids_above)
                            push!(watersurface,minpoint_high[2]) # store it as a watersurface point
                        else
                            continue
                        end

                    else    #if it is BELOW the mean
                        ids_below = ranges(rang[1],rang[2],voxel)
                        if length(ids_below) > 0
                            minpoint_low = findmin(ids_below) #1. height, 2. index in the voxel
                            inten = coords_dict[minpoint_low[2]][4] #intensity
                            imin[minpoint_low[2]] = [minpoint_low[1], counts[val], distance, inten]
                            else
                            continue
                        end
                    end

                elseif length(peaks) > 1

                    println("I have more than 1 peak")
                    println("length of bins ",length(bins))
                    println(counts, bins)

                    values = sort(counts[peaks], rev=true) #access the values of peaks, sort from the high to the low
                    print(values)

                    #1.
                    first = values[1]
                    getind_first = findall(x->x==first, counts)

                    #if the first HIGHEST peak is below the mean, then, it's the minpoint
                    if bins[getind_first[1]] <=  mean(bins)

                        density = first/npoints #density
                        distance = -1.0 #distance

                        rang_first = bins[getind_first[1]:getind_first[1]+1] #range in z level o f a bin
                        ids_fir = ranges(rang_first[1],rang_first[2],voxel)

                        if length(ids_fir) > 0
                            minpoint_low = findmin(ids_fir)
                            inten = coords_dict[minpoint_low[2]][4]
                            imin[minpoint_low[2]] = [minpoint_low[1], first, distance, inten]
                        else
                            continue
                        end
                    else #the first HIGHEST is above the mean, so watersurface!

                        #1. watersurface point
                        rang_first = bins[getind_first[1]:getind_first[1]+1]
                        ids_fir = ranges(rang_first[1],rang_first[2],voxel)

                        if length(ids_fir) > 0
                            minpoint_high = findmin(ids_fir)
                            push!(watersurface,minpoint_high[2])
                        else
                            continue
                        end
                        #2.
                        sec = values[2]
                        getind_sec = findall(x->x==sec, counts)

                        density = sec/npoints #density
                        distance = abs(bins[getind_first[1]] - bins[getind_sec[1]]) #distance

                        rang_sec = bins[getind_sec[1]:getind_sec[1] + 1]
                        ids_sec =  ranges(rang_sec[1],rang_sec[2],voxel)

                        if length(ids_sec) > 0
                            minpoint_low = findmin(ids_sec)
                            inten = coords_dict[minpoint_low[2]][4]
                            imin[minpoint_low[2]] = [minpoint_low[1], sec, distance, inten]
                        else
                            continue
                        end
                    end
                else
                    continue
                end
                println()
                #To SAVE the plots
                path = "C:/Users/alexandr/OneDrive - Stichting Deltares/Desktop/Thesis/P3_processing/with_greenstowa/Working_to/199NL3_Voxels/"
                path = joinpath(path,"nbins_" * string(length(bins)-1) * "_voxel_" * string(i))
                # plt.savefig(path,dpi=300)
                plt.close()
                global c+=1
            else
                imin[i] = [0, 0, 0, 0, 0] #if the number of bins is lower than 2, then I initialize all zeros
            end
        end
    end
end

#Get MAXIMUM values of density, distance, intensity for normalization procedure
den=[]; dis=[]; int=[];
for id in collect(values(imin))
    push!(den,id[2]); push!(dis,id[3]); push!(int,id[4])
end
maxden = maximum(den); maxdis = maximum(dis); maxint = maximum(int)

#Compute the Confidence values
confiden_vals = Dict{Int64,Array{Float64,1}}()
for i in keys(imin)
    if imin[i][3] == -1
        confiden_vals[i] = [imin[i][2]/maxden, 0, imin[i][4]/maxint]
    else
        confiden_vals[i] = [imin[i][2]/maxden, imin[i][3]/maxdis, imin[i][4]/maxint]
    end
end

#Define the threshold values of density, intensity and distance
nden=[]; ndis=[]; nint=[];
for id in collect(values(confiden_vals))
    push!(nden,id[1]); push!(ndis,id[2]); push!(nint,id[3])
end
#Mean and Median thresholds for the density, intensity and distance
mean_nden = mean(nden); mean_ndis = mean(ndis) ; mean_nint = mean(nint)
median_nden = median(nden); median_dis = median(ndis); median_nint = median(nint)

#write confidence values PER point in a las file
writeoutput = joinpath(workdir, mydict["dataname"] * mydict["outputname"] * ".laz")
write_confidence(writeoutput, ds, confiden_vals, median_nden, median_nint, mean_ndis)



## ------Write a RASTER --------- ##
np = pyimport("numpy")
rasterio = pyimport("rasterio")

#Create the extend of the array
keyes = collect(keys(vox_dict))
xxtend = Int64(round((nb_max[1] - nb_min[1])/voxelsize[1]))
yxtend = Int64(round((nb_max[2] - nb_min[2])/voxelsize[1]))
bands = 5 #number of bands
arr = np.zeros((xxtend,yxtend,bands)) #array in the bound box extend

#get confidence values per point
myconf = confidence(ds, confiden_vals, median_nden, median_nint, mean_ndis)

#store data into a multidimensional array
all = []
for key in keyes
    if key in collect(keys(imin))
        println("key: ",key)
        row, col, dontneedit = vox_dict[key]
        if row!=0 && col!=0 && dontneedit!=0
            arr[row,col,1] = imin[key][1] #z
            arr[row,col,2] = imin[key][2] #density
            arr[row,col,3] = imin[key][3] #distance
            arr[row,col,4] = imin[key][4] #intensity
            arr[row,col,5] = myconf[key][5] #confidence
            push!(all,(row,col))
        else
            continue
        end
    end
end


#Write raster output with the multidimensional array
rot = np.rot90(arr, k=1, axes=(0, 1)) #rotate the matrix
outtiff = joinpath(workdir,mydict["dataname"] * mydict["rasteroutput"])
f = rasterio.open(outtiff,
                  "w+",
                  driver= "GTiff",
                  width= xxtend,
                  height= yxtend,
                  count= bands,
                  dtype= rasterio.float64,
                  crs = "EPSG:28992",
                  nodata= -9999,)
f.write(rot[:,:,1],1)
f.write(rot[:,:,2],2)
f.write(rot[:,:,3],3)
f.write(rot[:,:,4],4)
f.write(rot[:,:,5],5)
f.close()
