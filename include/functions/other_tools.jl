function strToFile(fileName, string, mode) 
    # Appends "string" to the end of a file "fileName" and immediately closes it.
    # fileName and string must be strings
    # mode: r, read; w, write; a, append; r+, special read & write 
    open(fileName, mode) do f
        write(f, string)
    end
end
