
const myLightColors = [
  colorant"#FF74D4",
  colorant"#5386E4",
  colorant"#7BE0AD",
  colorant"#4B174B",
  colorant"#EE6C4D",
  colorant"#DBFE87",
  colorant"#716A5C",
  colorant"#D0A819"
]

const myLightTheme =  Theme(
  fontsize = 16,

  fonts = (
    regular     = "Verdana",
    bold        = "Verdana Bold",
    italic      = "Verdana Italic",
    bold_italic = "Verdana Bold Italic"
  ),

  Axis = (
    xlabelsize = 20,
    ylabelsize = 20,
    xgridstyle = :dash,
    ygridstyle = :dash
  ),

  Legend = (
    labelsize = 18,
  ),

  Lines = (
    linewidth = 2,
  ),

  palette = (patchcolor=collect(myLightColors), color=myLightColors)
)


function myAxis(f)
  ax = Axis(f)

  ax.xlabelsize = 18
  ax.ylabelsize = 18
  ax.titlesize  = 20

  ax.xlabelfont = "Verdana"
  ax.ylabelfont = "Verdana"
  ax.titlefont  = "Verdana Bold"

  return ax
end

# Theme(
#         resolution=RESOLUTION .* scale,
#         figure_padding=1,
#         rowgap=0,
#         colgap=0,
#         fontsize=fontsize,
#         fonts=(;regular="MinionPro-Capt"),
#         Axis=(
#             xminorticksvisible=true,
#             yminorticksvisible=true,
#             xgridvisible=false,
#             ygridvisible=false,
#             xtickalign=1,
#             ytickalign=1,
#             xminortickalign=1,
#             yminortickalign=1,
#             spinewidth=scale,
#             xtickwidth=scale,
#             xminortickwidth=scale,
#             ytickwidth=scale,
#             yminortickwidth=scale,
#             xminorticksize=2scale,
#             xticksize=3scale,
#             yminorticksize=2scale,
#             yticksize=3scale,
#             xlabelpadding=scale,
#             ylabelpadding=scale,
#         ),
#         palette=(patchcolor=collect(color), color=color),
#         Lines=(
#             linewidth=1.5scale,
#         ),
#         Scatter=(
#             strokewidth=0.75scale,
#             markersize=markersize*scale
#         ),
#         Legend=(
#             framevisible=true,
#             colgap=5scale,
#             rowgap=0,
#             patchsize=(10scale, 10scale),
#             patchlabelgap=3scale,
#             padding=(3scale, 3scale, 2scale, 2scale),
#             merge=true,
#             labelsize=fontsize,
#             linewidth=1.5scale,
#             markersize=markersize*scale,
#             markerstrokewidth=0.75scale,
#             titlesize=8scale,
#             titlegap=0
#         ),
#         Colorbar=(
#             tickalign=scale,
#             ticksize=3scale,
#             lip_vertical_label=true
#         ),
#     )