
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