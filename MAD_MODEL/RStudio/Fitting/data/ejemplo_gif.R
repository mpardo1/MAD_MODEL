# instalar librerias si es preciso
install.packages('gganimate') 
install.packages('ggplot2')
library(gganimate)
library(ggplot2)


set.seed(123)
data <- data.frame(x = rnorm(10000), y = rnorm(10000), t = 1:10000) # generar datos

# jerga ggplot, no se si estas familiarizada
giff <- ggplot(data = data, aes(x, y)) + # x, y van a ser posiciones
        geom_point() + # los datos van a ser representados en forma de puntos
        transition_time(t) +# ahi es donde va la variable tiempo, que hara variar las posiciones x y
        labs(title = 'Tiempo :{frame_time}') # para trazar el tiempo en el giff

animate(giff, duration = 10, fps = 30) # la duracion en segundos y los fps son parametros q puedes controlar

## renderizando ...

# ahi esta

# Para grabar el gif, le das right click en el "viewer", donde esta el gif, y le das a "Save Image"
# Ya veras que la extension del archivo es .gif por defecto.

