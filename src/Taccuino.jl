module Taccuino

using DataFrames
using DataArrays
using MAT
using CSV
using Dates

export get_data, createfilelist, paths_dataframe, convert2Bool, convert2Int, get_sessionname,
    get_CAMmousedate, get_BHVmousedate, get_Protocollo, get_streak, get_streakstart, get_sequence,
    get_correct, preprocess, create_pokes_single_session, create_pokes_dataframe, create_streak_dataframe,
    check_fiberlocation


"""
`createfilelist`
use get_data function to obtain all filenames of behaviour
"""

    function createfilelist(Directory_path::String, Mice_suffix::String)
        bhv = get_data(Directory_path,:bhv)
        #use get_sessionname to select relevant session (for instance use exp naming code)
        bhv_session = map(t->get_sessionname(t,Mice_suffix),bhv)#to be changed for each dataframe
        # get_sessionname return a start result for sessions that don't match the criteria this can be used to prune irrelevant paths
        bhv = bhv[bhv_session.!="start"]
        bhv_session = bhv_session[bhv_session.!="start"]
        return bhv
    end

"""
`get_data`

Functions to find the path to data files according to character pattern considering that each session is a subfolder
Get data revised: this version operate in 2 possible way,
in the first way it will looks for file names containing a specified string
in the second way use a symbol to refer to a dictionary and find the specified string
"""
    #Method 1 inputs a directory and a string
    function get_data(dirnames,what::String)
        location = String[] #array to load with all the paths corrisponding to researched file type
        if eltype(dirnames)== Char #in case only one folder is loaded the for loop would research in single character
            tool = String[]
            push!(tool,dirnames)
            dirnames = tool
        end
        for dirname in dirnames
            files = readdir(dirname)
            for file in files
                if ismatch(Regex(what), file)
                    complete_filename = joinpath(dirname,file)
                    push!(location,complete_filename)
                end
            end
        end
        return location
    end
    #Method 2 a directory and a symbol
    function get_data(dirnames, kind::Symbol)
        #the dictionary refers the symbol in the input to a specific string to look for
        ext_dict = Dict(:bhv => "a.csv", :cam => ".mat", :log => "AI.csv")
        if !(kind in keys(ext_dict))
            error("Choose among $(keys(ext_dict))")
        end
        #once the string is identified the function call itself again with the first method
        return get_data(dirnames, ext_dict[kind])
    end

"""
`get_sessionname(filepath, what::String)`

Use it to find the name of a session from a path, can accept a string or a symbol connected to a dict to find
file matching the requirements
"""
#this function extract the name of a session from a filepath according to a given character pattern
function get_sessionname(filepath, what::String)
    pathinfo = split(filepath,"/")
    sessionname = "start"
        for piece in pathinfo
            if ismatch(Regex(what), string(piece))
                sessionname = piece
            end
        end
    return sessionname
end

# the second method allows to save experiment related character pattern in a dictionary
function get_sessionname(filepath, kind::Symbol)
    #the dictionary refers the symbol in the input to a specific string to look for
    ext_dict = Dict(:GcAMP => "170", :BilNac => "NB")
    if !(kind in keys(ext_dict))
        error("Choose among $(keys(ext_dict))")
    end
    #once the string is identified the function call itself again with the first method
    return get_sessionname(filepath, ext_dict[kind])
end

"""
`paths_dataframe`
Create a Dataframe to store paths of files to preprocess
"""
function paths_dataframe(bhv)
    behavior = DataFrame()
    behavior[:Bhv_Path]= bhv
    ##### extract date and mouse ID per session using get_mousedate (it works with a full path)
    MouseID = Array{String}(size(behavior,1))
    Day = Array{String}(size(behavior,1))
    Session = Array{String}(size(behavior,1))
    for i = collect(1:size(behavior,1))
        MouseID[i], Day[i], Session[i] = get_BHVmousedate(behavior[i,:Bhv_Path])
    end
    behavior[:MouseID] = MouseID
    behavior[:Day] = Day#file properties are not reliable for the date of the session
    behavior[:Session] = Session.*".csv";
    return behavior
end


"""
`convert2Bool(df,symbols)`

Converts values of columns symbols of dataframe df in Bool
"""
function convert2Bool(df,symbols)
    for symbol in symbols
        df[:,symbol] = Bool.(df[:,symbol])
    end
end

"""
`convert2Int(df,symbols)`

Converts values of columns symbols of dataframe df in Int
"""
function convert2Int(df,symbols)
    for symbol in symbols
        df[:,symbol]=Int64.(df[:,symbol])
    end
end

"""
`get_CAMmousedate(filepath, pattern)`

it extract the session name by pattern match:
It assumes that the date of the creation of the file is the day of the session
"""
function get_CAMmousedate(filepath, pattern)
    fileinfo = get_sessionname(filepath,pattern)
    mouse, note = split(fileinfo, "_")[1:2]# decompose the file name by _ and take the first as animal name
    #and the second as extra experimental note, like target area
    date = Dates.Date(Dates.unix2datetime(ctime(filepath)))#return the date from file properties,
    #convert it from unix to normal time and take only the date
    return mouse, date, note
end

"""
`get_BHVmousedate(filepath)`

it extract the session name by Regular expression match:
It assumes that session name is composed by
2 characters and a number for the mouse
followed by the date
"""
function get_BHVmousedate(filepath)
    sessionREGEX = match(r"[a-zA-Z]{2}\d+_\d{6}",filepath); #the result is a regex object with several info
    session = sessionREGEX.match;
    mouse, giorno = split(session,"_")[1:2]
    day = "20"*giorno
    return mouse, day, session
end
"""
`get_Protocollo(df)`

Python preprocessing creates a Boolean vector to distinguish between 2 Protocollo
relative Gamma and Probability are infered by the respective GammaVec and ProbVec
"""
function get_Protocollo(df)
    ProtName = String[]
    for i = collect(1:size(df,1))
        #in the original data protocol are labeled either 0 or 1
        #this function collapse the respective gamma and prob in a string
        if df[i,:Protocollo]==0
            curr_prot = string(df[i,:ProbVec0],"/",df[i,:GamVec0])
            push!(ProtName,curr_prot)
        elseif df[i,:Protocollo]==1
            curr_prot = string(df[i,:ProbVec1],"/",df[i,:GamVec1])
            push!(ProtName,curr_prot)
        end
    end
    df[:Protocol] = pool(ProtName)
end
"""
`get_streak(df)`

Starting from 1 create a counter that increase
when detect change in side between a poke and the previous
"""
function get_streak(df)
    Streak = Int64[1] #create an array to fill with streak counter first value is by definition streak 1
    for i = 2:size(df,1)
        if df[i,:Side] != df[i-1,:Side] #if previous side is different from current side
            push!(Streak, Streak[i-1]+1) #increase the counter
        else
            push!(Streak, Streak[i-1]) #otherwise keep the counter fixed
        end
    end
    return Streak
end
"""
`get_sequence`

Starting from 1 create a counter that increase when detect change in  a categorical variable
a 3rd argument can be use to reset the counter at the change of another categorical variable
"""
function get_sequence(df,category,by)
    Streak = Int64[1] #create an array to fill with streak counter first value is by definition streak 1
    for i = 2:size(df,1)
        if df[i-1,by] == df[i,by]
            if df[i,category] != df[i-1,category] #if previous side is different from current side
                push!(Streak, Streak[i-1]+1) #increase the counter
            else
                push!(Streak, Streak[i-1]) #otherwise keep the counter fixed
            end
        else
            push!(Streak, 1)
        end
    end
    return Streak
end
##Method2 giving dataframe
function get_sequence(df,category::Symbol)
    Streak = Int64[1] #create an array to fill with streak counter first value is by definition streak 1
    for i = 2:size(df,1)
        if df[i,category] != df[i-1,category] #if previous side is different from current side
            push!(Streak, Streak[i-1]+1) #increase the counter
        else
            push!(Streak, Streak[i-1]) #otherwise keep the counter fixed
        end
    end
    return Streak
end


"""
`get_streakstart(df)`

function to create a boolean column that flags the begin of a streak
begin of trials
"""
function get_streakstart(df)
    #create an array to point first poke of a given trial,
    #by definition first poke of the session is the beginning of a trial
    Streakstart = Bool[true]
    for i= 2:size(df,1)
        #if previous streak counter is different from the current
        if df[i,:StreakCount] != df[i-1,:StreakCount]
            push!(Streakstart,true) #begin of a new trial
        else
            push!(Streakstart,false)#otherwise same trial
        end
    end
    return Streakstart
end

"""
`get_correct`

TO BE REDEFINE WORKS ONLY FOR 100% GAMMA NOT SO WELL also
function to define correct and incorrect trials
"""
function get_correct(df) #this function works correctly only for protocols with 100% probability flipping gamma
    Correct = Bool[] #create an array to fill with streak counter
    for i = 1:size(df,1)
        if df[i,:Side] == df[i,:SideHigh] #if current side is equal from side high
            push!(Correct, true) # set correct true
        else
            push!(Correct, false) #otherwise false
        end
    end
    return Correct
end

"""
`preprocess`

function to preprocess flipping behavioural data
from the python csv file combining the previous function
"""

function preprocess(bhv_files)
    curr_data=CSV.read(bhv_files,nullable = false)
    rename!(curr_data, Symbol(""), :PokeCount) #change poke counter name
    curr_data[:PokeCount]= curr_data[:PokeCount]+1
    booleans=[:Reward,:Side,:SideHigh,:Stim] #columns to convert to Bool
    integers=[:Protocollo,:ProbVec0,:ProbVec1,:GamVec0,:GamVec1,:delta] #columns to convert to Int64
    convert2Bool(curr_data,booleans)
    convert2Int(curr_data,integers)
    mouse, day, session = get_BHVmousedate(bhv_files)
    curr_data[:PokeDur] = curr_data[:PokeOut]-curr_data[:PokeIn]
    ##### filter invalid pokes
    curr_data[:MouseID] = mouse
    curr_data[:Day] = day
    curr_data[:Session] = session
    get_Protocollo(curr_data)#create a columns with a unique string to distinguish protocols
    curr_data[:StreakCount] = get_sequence(curr_data,:Side)
    curr_data[:StreakStart] = get_streakstart(curr_data)
    curr_data[:Correct] = get_correct(curr_data)
    for x in[:ProbVec0,:ProbVec1,:GamVec0,:GamVec1,:Protocollo]
        delete!(curr_data, x)
    end
    return curr_data, session
end

"""
`check_fiberlocation`

look for a dataset where fiberlocation across day is stored
"""
function check_fiberlocation(data,Exp_name)
    filetofind=joinpath("/Users/dariosarra/Google Drive/Flipping/run_task_photo/"*Exp_name*"/FiberLocation"*".csv");
    if isfile(filetofind)
        fiberlocation = readtable(filetofind);
        pokes_table = join(data, fiberlocation, on = :Session, kind = :inner);
        println("found fibres location file, HAVE YOU UPDATED IT?")
    else
        println("no fibres location file")
        pokes_table = data;
    end
    return pokes_table
end

"""
`create_pokes_single_session`

Preprocess all the selected files and create single file per session
"""
##### Preprocess all the selected files and create single file per session
function create_pokes_single_session(behavior::DataFrame, Exp_name::String)
    c=0
    b=0
    Preprocessed_path = []; #to push a string is better to create an empty vector
        for i=1:size(behavior,1)
        path = behavior[i,:Bhv_Path]
        session = behavior[i,:Session]
        filetosave=joinpath("/Users/dariosarra/Google Drive/Flipping/run_task_photo/"*Exp_name*"/Bhv/"
                *session)
        push!(Preprocessed_path,filetosave)
            if ~isfile(filetosave)
                data,bho = preprocess(path)
                writetable(filetosave,data)
                b=b+1
            else
                c=c+1
            end
        end
    println("Existing file = ",c," Preprocessed = ",b)
    behavior[:Preprocessed_Path] = Preprocessed_path;
    return behavior
end

"""
`create_pokes_dataframe`
join all the preprocessed pokes dataframe in a single dataframe
"""
function create_pokes_dataframe(behavior::DataFrame,Exp_type::String,Exp_name::String)
    case = 0;
    data = DataFrame()
    for i in behavior[:Preprocessed_Path]
        if case == 0;
            data = readtable(i);
            case = 1;
            elseif case ==1;
            x = readtable(i);
            append!(data,x);
        end
    end
    #### save the pokes dataframe
    filetosave = "/Users/dariosarra/Google Drive/Flipping/Datasets/"*Exp_type*"/"*Exp_name*"/pokes"*Exp_name*".csv"
    writetable(filetosave,data)
    return data
end

"""
`create_streak_dataframe`
From the pokes dataframe creates one for streaks
"""
function create_streak_dataframe(data::DataFrame,Exp_type::String,Exp_name::String)
    streak_table = by(data, [ :Day, :MouseID, :StreakCount,]) do df
        DataFrame(Wall = df[1,:Wall],
                  Side = df[1,:Side],
                  SideHigh = df[1,:SideHigh],
                  Stim = df[1,:Stim],
                  Num_pokes = size(df,1),
                  Num_Rewards = length(find(df[:Reward].==1)),
                  Last_Reward = findlast(df[:Reward] .==1),
                  Trial_duration = (df[:PokeOut][end]-df[:PokeIn][1]),
                  Start = (df[1,:PokeIn]),
                  Stop = (df[end,:PokeOut]),
                  Correct = df[1,:Correct],
                  Session = df[1,:MouseID]*"_"*string(df[:Day][1]),
                  Session2 = df[1,:Session],
                  Area = df[1,:Area],
                  Protocollo = df[1,:Protocol]
                  )
              end

    streak_table[:AfterLast] = streak_table[:Num_pokes] - streak_table[:Last_Reward];
    sort!(streak_table, cols = [order(:Session), order(:StreakCount)])


    #travel duration
    streak_table[:Travel_duration] = Array{Float64}(size(streak_table,1));
    a= streak_table[:Stop][2:end] - streak_table[:Start][1:end-1];
    streak_table[1:end-1,:Travel_duration] = (streak_table[:Stop][2:end] - streak_table[:Start][1:end-1]);
    streak_table[end,:Travel_duration] = 0;

    #block counter
    streak_table[:Block_counter] = get_sequence(streak_table,:Wall,:Session);#va diviso per sessions;

    #Block length
    streak_table[:Block_Streak] = get_sequence(streak_table,:StreakCount,:Block_counter);

    #### Save streaktable
    filetosave = "/Users/dariosarra/Google Drive/Flipping/Datasets/"*Exp_type*"/"*Exp_name*"/streaks"*Exp_name*".csv"
    writetable(filetosave,streak_table)
    return streak_table
end
end
