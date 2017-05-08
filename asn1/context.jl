
"
Describe the reading context.
Which element types are expected for a compound type (SEQUENCE, SET)?
Which ordinal number (SEQUENCE)?
Which element type is expected ((SEQUENCE, SET) OF)?
Which tags are included?
"
abstract Expected

immutable Context
	tag::Tag
	utc::Expected
end

immutable ExpectedNull <: Expected
end

const EXPECTED_NULL = ExpectedNull()
const CONTEXT_NULL = Context(Tag(UT_EOC), EXPECTED_NULL)

immutable ExpectedSequence <: Expected
	elementType::Vector{Context}
	ordinal::Int	# position in sequence definition
end

immutable ExpectedSequenceOf <: Expected
	elementType::Context
end

immutable ExpectedSet <: Expected
	elementTypes::Vector{Context}
end

immutable ExpectedSetOf <: Expected
	elementType::Context
end

function subcontext(utc::ExpectedSequence)::Context
	utc.ordinal <= length(utc.elementTypes) ? utc.elementTypes[utc.ordinal] : CONTEXT_NULL
end

function nextcontext(utc::ExpectedSequence)::Expected
	utc.ordinal < length(utc.elementTypes) ? ExpectedSequence(utc.elementTypes, utc+1) : EXPECTED_NULL
end

subcontext(utc::ExpectedSet)::Context = utc
subcontext(utc::ExpectedSequenceOf)::Context = utc.elementType 
subcontext(utc::ExpectedSetOf)::Context = utc.elementType 
subcontext(utc::Expected)::Context = CONTEXT_NULL
nextcontext(utc::Expected)::Expected = utc

subcontext(ctx::Context)::Context = subcontext(ctx.utc)
function nextcontext(ctx::Context)::Context
	utc = nextcontext(ctx.utc)
	utc == ctx.utc ? ctx : utc == EXPECTED_NULL ? CONTEXT_NULL : Context(ctx.tag, utc)
end

"tail function returned if no matching tail function exists"
function errortail(tag::Tag, ctx::Context)
	failparse("no tail function for tag $tag")
end

"derive type of data following in the input stream from the peceding tag
and the reading context"
function tailfunction(tag::Tag, ctx::Context)::Function
	iseoc(ctx.tag) || ctx.tag == tag || failparse("expected tag $(ctx.tag) - got $tag")
	get(TAILS, tag, errortail)
end
