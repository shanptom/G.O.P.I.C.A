# MetaPiX Shiny Application
# This is the main entry point for the MetaPiX application.
# It sources the UI and server logic from separate files and launches the app.

# Source the UI and Server components
source("ui.R")
source("server.R")

# Launch the Shiny App
shinyApp(ui = ui, server = server)
