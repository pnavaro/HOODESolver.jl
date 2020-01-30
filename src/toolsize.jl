
# Helper: Pointer cache is used to break circular references in objects and avoid double countings
is_seen!(x, ptr_cache) = in(pointer_from_objref(x), ptr_cache) || (push!(ptr_cache, pointer_from_objref(x)); false)

# Helper: Catch types without size method
sizeof_catch(x) = try sizeof(x) catch nothing end

# Composite Types
function totalsizeof(x, ptr_cache = Set())
#  is_seen!(x, ptr_cache) && return 0
  result = sizeof_catch(x)

  this_names = nothing

  # Check if type has names method
  try
    this_names = names(x)
  catch
    this_names = []
  end

  el = nothing
  isempty(this_names) ? result :
    result + mapreduce(i -> 
      # Check if some elements are not reachable
      (try el = getfield(x,i) catch nothing end; totalsizeof(el, ptr_cache))
      , +, this_names)
end

# Tuples
function totalsizeof(x::Tuple, ptr_cache = Set())
  # is_seen!(x, ptr_cache) && return 0
  isempty(x) ? 0 : mapreduce(i -> totalsizeof(i, ptr_cache), +, x)
end

# Array
function  totalsizeof(x::Array{T}, ptr_cache = Set()) where{T<:Any}
 # is_seen!(x, ptr_cache) && return 0
  result = sizeof(x)
  el = nothing
  for i in 1:length(x)
    # Check if some elements are undefined
    try el = x[i] catch  nothing end
    result += totalsizeof(el, ptr_cache)
  end
  result
end

# String, Function, IntrinsicFunction
totalsizeof(x::Union{String, Function}, ptr_cache = Set()) = is_seen!(x, ptr_cache) ? 0 : sizeof(x)
