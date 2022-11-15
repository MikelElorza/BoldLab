### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 80a4e1cc-4d59-4c8a-9d65-1fe0d2d77bf2
using Pkg; Pkg.activate(ENV["JBoldLab"])

# ╔═╡ 23b4e644-8331-4aa6-8369-47ce3ff0e143
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Images
	using ImageBinarization
	using Colors
	using Plots
	using Printf
	using Interpolations
	using QuadGK
	using Markdown
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
	using Peaks
	using FFTW
	using DSP
	using Clustering
	import Glob
end

# ╔═╡ 3a7182c2-8ffd-4d6a-975f-42f9ff1b0744
using Test

# ╔═╡ 9f71bc31-df8f-48c5-815a-c8db9e46100f
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 19035c0d-93a7-4bb9-be05-d2a9b9ac4619
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ f4e379c4-a2b2-4703-bcbc-f6d7c996354a
lbl = ingredients("../src/BoldLab.jl")

# ╔═╡ d88ba94f-dd9e-4c46-adac-4c3f7e556cad
lbl.BoldLab.print_greeting()

# ╔═╡ 51710bb8-9c4a-4abc-9fe9-02e87bd4e4c5
PlutoUI.TableOfContents(title="Laser Lab CMOS analysis", indent=true)

# ╔═╡ d1ace15e-fe1a-4552-9144-a0824ae8ae0f
md"""
# Contributors

- J.J. Gómez Cadenas (jjgomezcadenas@gmail.com)
- M. Elorza
- P. Herrero
"""

# ╔═╡ 0122a81d-ee51-4355-8ddb-6429f15e8ea1
md"""
# Code
"""

# ╔═╡ dcd1a13d-278c-4168-9649-c69ccfe34633
md"""
## Dropbox path
Here you must enter de path of your LaserLab Dropbox. With this path, the program can access to every file in the Dropbox.
"""

# ╔═╡ 10cf446d-4fa9-40c3-99c7-9160f6bcd25c
LaserLabp="C:\\Users\\Mikel\\LaserLab Dropbox\\"

# ╔═╡ e2797851-5f39-4f82-b98c-916975f4ed02
md""" 
## Select run
"""

# ╔═╡ a0842787-f875-46f0-a47c-af5e792110c9
begin
bodir=joinpath(LaserLabp,"Proyectos\\pdata\\")
bcmdir=joinpath(LaserLabp,"Proyectos\\data\\")
end

# ╔═╡ 2f5d32e7-1196-4ee3-83fd-9ccd23467639
begin
cmds = filter(name->occursin("CMOS",name), readdir(bcmdir))
md""" Select CMOS dir : $(@bind scmos Select(cmds))"""
end

# ╔═╡ b686f66c-ab43-4185-ace0-88d8e4516c45
cmdir=joinpath(bcmdir,scmos)

# ╔═╡ 853b8b2e-66ed-4723-8884-213e5fd4a0e7
md"""
# Tests
"""

# ╔═╡ 20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
md"""
# Functions
"""

# ╔═╡ ff60879d-c013-4b4e-b121-e6cfa3f5517c
function select(dir)
	readdir(dir)
end

# ╔═╡ 74e548a7-4d15-46a8-9819-95f411910115
let
subs = select(cmdir)
md""" Select substrate : $(@bind ssub Select(subs))"""
end

# ╔═╡ a9e40c19-42e6-4a4d-ab31-c4ff5ae3d75a
subsp=joinpath(cmdir,ssub)

# ╔═╡ cf00e1e9-35d3-480b-a8bb-f842db0ac7b7
let
exps = select(subsp)
md""" Select substrate : $(@bind sexp Select(exps))"""
end

# ╔═╡ fa13aaf4-6894-4360-bb28-e425f9392316
expp=joinpath(subsp,sexp)

# ╔═╡ 0596c3df-2bba-4324-b99b-494fbdf779c6
let
runs = select(expp)
md""" Select run : $(@bind srun Select(runs))"""
end

# ╔═╡ a7e1bb4d-1e92-419f-b609-43d5384a6685
runp=joinpath(expp,srun)

# ╔═╡ e14f380d-9f64-4efd-9b40-07e7105d8704
@test typeof(select(runp))==Vector{String}

# ╔═╡ 0bf5c4a7-b369-4438-a987-253f242dd88f
function dummy(x,y)
	x+y
end

# ╔═╡ 71ed6fbd-0342-4cf5-9d8c-aa8f791d85f1
@test dummy(3,4) == 7

# ╔═╡ Cell order:
# ╠═80a4e1cc-4d59-4c8a-9d65-1fe0d2d77bf2
# ╠═23b4e644-8331-4aa6-8369-47ce3ff0e143
# ╠═3a7182c2-8ffd-4d6a-975f-42f9ff1b0744
# ╠═9f71bc31-df8f-48c5-815a-c8db9e46100f
# ╠═19035c0d-93a7-4bb9-be05-d2a9b9ac4619
# ╠═f4e379c4-a2b2-4703-bcbc-f6d7c996354a
# ╠═d88ba94f-dd9e-4c46-adac-4c3f7e556cad
# ╠═51710bb8-9c4a-4abc-9fe9-02e87bd4e4c5
# ╠═d1ace15e-fe1a-4552-9144-a0824ae8ae0f
# ╠═0122a81d-ee51-4355-8ddb-6429f15e8ea1
# ╠═dcd1a13d-278c-4168-9649-c69ccfe34633
# ╠═10cf446d-4fa9-40c3-99c7-9160f6bcd25c
# ╠═e2797851-5f39-4f82-b98c-916975f4ed02
# ╠═a0842787-f875-46f0-a47c-af5e792110c9
# ╠═2f5d32e7-1196-4ee3-83fd-9ccd23467639
# ╠═b686f66c-ab43-4185-ace0-88d8e4516c45
# ╠═74e548a7-4d15-46a8-9819-95f411910115
# ╠═a9e40c19-42e6-4a4d-ab31-c4ff5ae3d75a
# ╠═cf00e1e9-35d3-480b-a8bb-f842db0ac7b7
# ╠═fa13aaf4-6894-4360-bb28-e425f9392316
# ╠═0596c3df-2bba-4324-b99b-494fbdf779c6
# ╠═a7e1bb4d-1e92-419f-b609-43d5384a6685
# ╠═853b8b2e-66ed-4723-8884-213e5fd4a0e7
# ╠═71ed6fbd-0342-4cf5-9d8c-aa8f791d85f1
# ╠═e14f380d-9f64-4efd-9b40-07e7105d8704
# ╠═20770d2f-ca8f-4fb3-b11d-d00f93e3a0cc
# ╠═ff60879d-c013-4b4e-b121-e6cfa3f5517c
# ╠═0bf5c4a7-b369-4438-a987-253f242dd88f
