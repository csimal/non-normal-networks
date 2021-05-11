
function si!(du,u,p,t)
    β = p
    s = @view u[:,1]
    i = @view u[:,2]
    @. du[:,1] = -β * s * i
    @. du[:,2] = - du[:,1]
end

function sis!(du,u,p,t)
    β, γ = p[1], p[2]
    s = @view u[:,1]
    i = @view u[:,2]
    @. du[:,1] = (-β * s * i) + (γ * i)
    @. du[:,2] = - du[:,1]
end

function sir!(du,u,p,t)
    β, γ = p[1], p[2]
    s = @view u[:,1]
    i = @view u[:,2]
    @. du[:,1] = -β * s * i
    @. du[:,2] = - du[:,1] -  (γ * i)
    @. du[:,3] = γ * i
end

function ()
    
end