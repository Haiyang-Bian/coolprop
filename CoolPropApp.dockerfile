FROM julia
ENV MYPATH = /home/coolprop
WORKDIR ${MYPATH}
RUN julia -e "using Pkg; Pkg.add([\"CoolProp\",\"Unitful\",\"Oxygen\"])"
COPY . .
EXPOSE 8080
VOLUME [ "/data" ]
CMD [ "julia", "Router.jl" ]