def decorator_function(original_function):
    def wrapper_function():
        print('wrapper executed before {}'.format(original_function.__name__)) # do this thing
        return original_function() # return the output of the original function
    return wrapper_function # returns the wrapper function. so basically when you call the decorated function,
                            # it runs the original function

@decorator_function
def display():
    print('display function ran')

### above, this is exactly analogous to doing

display = decorator_function(display)