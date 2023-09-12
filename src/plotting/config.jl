
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
  fontsize = 18,

  fonts = (
    regular     = Makie.texfont(),
    bold        = Makie.texfont(:bold),
    italic      = Makie.texfont(:italic),
    bold_italic = Makie.texfont(:bolditalic)
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
