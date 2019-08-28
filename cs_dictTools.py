def populateDictOfLists(mydict, queryKey, queryValue):
    # check if key is already in dictionary
    if queryKey in mydict:
        # check if queryValue in list
        if queryValue in mydict[queryKey]:
            pass
        else:
            mydict[queryKey].append(queryValue)

    else:   # if not, add it to the dictionary and run func again
        mydict[queryKey] = []
        populateDictOfLists(mydict, queryKey, queryValue)

    return mydict
