function include_lib( folder::String )
    for (root, _, files) in walkdir(folder)
        for file in files
            if endswith(file, ".jl")
                full_path = joinpath(root, file)
                # @info "Including $full_path"
                include(full_path)
            end
        end
    end
end
