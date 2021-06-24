using Documenter
using julia_app

makedocs(
    sitename = "julia_app",
    format = Documenter.HTML(),
    modules = [julia_app]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
