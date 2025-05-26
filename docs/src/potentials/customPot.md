# Custom Potentials

Using a custom potential within JMD is fairly straightforward. A custom potential with full functionality requires you make 5 functions, however, depending on your usecase you can make fewer. 

### Initializer Function

This function needs to initialize the variables used within your custom potential. The function should take a single parameter, either a vector of Atoms or a JMD cell object, and return a struct containing the potential's variables.
