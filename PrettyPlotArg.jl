
function PlotBasic()
    arg  = Dict(
         :framestyle => :box,
         :grid => false)

    return arg
end

function PrettyLeg(;fontsize::Int64=14,location::Symbol=:best)
    arg  = Dict(
        :legendfontsize => fontsize-2,
        :legendtitlefontsize => fontsize-2,
        :legend => location,
        :foreground_color_legend => :silver,
        :background_color_legend => :white
        )
    return arg
end

function PlotFontSize(;fontsize::Int64=14)
    arg  = Dict(
        :titlefontsize => fontsize, 
        :xguidefontsize => fontsize, 
        :yguidefontsize => fontsize, 
        :tickfontsize => fontsize-3)
    return arg 
end

