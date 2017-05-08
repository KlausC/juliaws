#! /bin/env julia

module MatlabTranslation

export translate

using ParserCombinator

function main(args...)
	argc = length(args)
    if argc == 0
        out = translate(STDIN)
    	for token in out
            @printf(ou, "%s", token)
        end
    end
    for fname in args
        in = open(fname, "r")
        ou = open(string(fname,".jl"), "a")
        try
    	    out = translate(in)	
    	    for token in out
                @printf(ou, "%s", token)
            end
        catch x
            println(x)
        finally
            close(in)
            close(ou)
        end
    end
end

function translate(input::IO)
    parse_lines(input, file)
end

# Transformers
function convertcomment(x)
    string('#', x, '\n')
end

function convertquoted(x)
    string('"', x, '"')
end

function convertwords(x)
    if ( x == "fprintf" )
        x = "@printf"
    end
    x
end

function pushargs(x...)
    push!(oargs, "")
    string(x...)
end

function addoptlist(x...)
    if length(x) == 1
        string(x[1], "()")
    else
        string(x...)
    end
end

function addend(x...)
    oa = pop!(oargs)
    if length(x) == 0 || ! endswith(x[1], "end")
        string(x..., oa, "\nend")
    else
        string(x..., oa)
    end
end

oargs = AbstractString[] 

function memoize_oargs(x...)
    if length(x) > 0
        oargs[end] = rstrip(lstrip(string(x...), '['), ']')
    end
    "" 
end

function addoutput(x...)
	if length(x) < 2
        string(x..., oargs[end])
    else
        string(x...)
    end
end


# Transformer Generators
function st1(name)
    function str(args...)
        string("<", name, ":", args..., ":", name, ">")
    end
    str
end

function st0(name)
    string
end

set_st!(s) = st = s
st = st0

# surround delimiters with optional spaces
del_start(p::Matcher) = Seq(p, spc)
del_middle(p::Matcher) = Seq(spc, p, spc)
del_end(p::Matcher) = Seq(spc, p)

# surround keywords delimiters with spaces
keyw_start(p::Matcher) = Seq(p, splus)
keyw_middle(p::Matcher) = Seq(splus, p, splus)
keyw_end(p::Matcher) = Seq(splus, p)

# Lexer
LF = p"(\r?\n)|\r"
comment = E"%" + p"[^\r\n]*" + Drop(LF) > convertcomment
lineend = LF | comment
contin  = Drop(E"..." + p"[ \t]*" + LF)
spp = Plus(p"[ \t]+" | contin)	# at least one space on same logical line
sp  = Opt(spp)	> string		# spaces on logical line
splus = Plus(lineend | spp) > string	# at least one space on same or more lines
spc = Opt(splus) > string		# spaces spread over one or more lines
PO = del_start(e"(")
PC = del_end(e")")
BO = del_start(e"[")
BC = del_end(e"]")
CO = del_start(e"{")
CC = del_end(e"}")
keyword = Alt(e"end", e"else", e"do", e"if", e"elseif", e"begin", e"function", e"return")

uop = p"[-+!~]"					# unary operators
number = p"(\.\d+|\d+\.?\d*)([eE][-+]?\d+)?i?"	# > st("number")
#number = rnumber + Opt(Opt(del_middle(p"[-+]") + rnumber) + e"i") > string
op = p"[.]?[-+\*/\\~?|><=&$§!^°]+"					# operator				
quoted = E"'" + p"[^']*" + E"'" > convertquoted
word = p"[[:alpha:]]\w*" > convertwords						# identifier
psep = del_middle(e",")	> string
bsep = spp | del_middle(p"[,;]" | lineend) > string
lsep = psep
sep = del_middle(e";" | lineend) > string


pop = p"\.?'" > string	# postscript unary operator
speci = e"." + word	> string									# fieldname specifier

expr = Delayed()
pstart = PO + expr
plist = pstart + psep + Star(expr + lsep) + Opt(expr) + PC	# tuple
flist = (PO + PC ) | (pstart + Star(psep + expr) + PC)	# function argument list
pexpr = pstart + PC										# parenthesized expression
blist = BO + expr + Star(bsep + expr) + Opt(bsep) + BC	# array value list
alist = BO + expr + Star(psep + expr) + BC				# array indexing list
olist = BO + Opt(word + Star(lsep + word)) + BC			# output argument list
llist = PO + Opt(word + Star(lsep + word)) + PC			# input argument list
dlist = flist | alist

prim = word | quoted | number | pexpr | plist | blist
term_tail = Star( pop | (sp + dlist) | (sp + speci) )
term = Star(uop + sp) + prim + term_tail
uexp = term + Star(del_middle(op) + term)	>st("uexp")
expr.matcher= (del_start(e":") + Opt(uexp)) | uexp +( del_middle(e":") +uexp)[0:2] > string

statement = Delayed()
statement_list = StarList(statement, Plus(sep)) + Star(sep)
oarglist = Alt(olist, word) + Drop(del_middle(e"=")) > memoize_oargs
funct_keyword = keyw_start(e"function") > pushargs
funct_headstart = Seq!(funct_keyword, Opt(oarglist), word) > string
funct_head = funct_headstart + Opt(spc + llist) > addoptlist
funct_end = Alt(keyw_end(e"end"), del_end(Eos())) > addend
funct = funct_head + spc + statement_list + funct_end > st("funct")

return_statement = keyw_start(e"return") + Opt(expr) > addoutput
block_statement = keyw_start(e"begin") + statement_list + keyw_end(e"end")
if_statement = Seq(keyw_start(e"if"), expr, splus, statement_list, 
               Star(keyw_middle(e"elseif") + statement_list),
			   Opt(keyw_middle(e"else") + statement_list), keyw_end(e"end"))
for_statement = keyw_start(e"for") + expr + Drop(Opt(psep)) + splus + statement_list + keyw_end(e"end") > st("for")
while_statement = keyw_start(e"while") + expr + splus + statement_list + keyw_end(e"end")

statement.matcher = Alt(
    return_statement,
    block_statement,
    if_statement,
    for_statement,
    while_statement,
    funct,
	keyword + Fail(),
	expr,
) > st("statement")

file = spc + Opt(statement_list) + Opt(sep) + Eos()

end #MatlabTranslation

using ParserCombinator
using MatlabTranslation
M = MatlabTranslation

#M.main()
