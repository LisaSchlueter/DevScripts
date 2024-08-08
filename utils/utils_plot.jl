
function _get_def_peakcolors()
    def_colors = Dict(:Tl208a => :dodgerblue, 
            :Bi212a => :darkorange, 
            :Tl208b => :lightseagreen,
            :Tl208DEP => :midnightblue,
            :Bi212FEP => :crimson,
            :Tl208SEP => :olive,
            :Tl208FEP => :violet)
    return def_colors
end

function _get_def_histargs(; fs = 14)
    HistArg = Dict(:normalize => :probability,
               :fill => :true,
               :foreground_color_legend => :silver,
               :background_color_legend => :white,
               :grid => :off,
               :xguidefontsize => fs+1, :xtickfontsize => fs-2,
               :yguidefontsize => fs+1, :ytickfontsize => fs-2,
               :legendfontsize => fs-2)
    return HistArg
end

function _def_plot(; fs::Integer = 14)
    default(foreground_color_legend = :silver,
    background_color_legend = :white,
    grid = :off,
    framestyle = :semi,
    xtickfontsize = fs-2,
    ytickfontsize = fs-2,
    legendfontsize = fs-2,
    xlabelfontsize = fs + 2,
    ylabelfontsize = fs + 2)
end
