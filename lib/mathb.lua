local round = function(number, digits)
    if (number == nil) then error("yeet",2) end
    digits = digits or 0
    number = number * (10^digits)
    number = number + 0.5
    number = math.floor(number)
    number = number / (10^digits)
    return number
end

local factorial = function( a ) 
    if (math.floor(a) ~= a or a < 0) then error("factorial: bad input",2) end
    if (a == 0) then return 1 end
    local output = 1
    for i=1,a do
        output = output * i
    end
    return output
end

return {
    round = round,
    factorial = factorial,
}
