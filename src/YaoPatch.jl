export Bag, add!

mutable struct Bag{N}<:TagBlock{AbstractBlock, N}
    content::AbstractBlock{N}
end

Yao.content(bag) = bag.content
Yao.chcontent(bag::Bag, content) = Bag(content)
Yao.mat(bag::Bag) = mat(bag.content)
Yao.apply!(reg::AbstractRegister, bag::Bag) = apply!(reg, bag.content)
YaoBlocks.PreserveStyle(::Bag) = YaoBlocks.PreserveAll()
setcontent!(bag::Bag, content) = (bag.content = content; bag)

function YaoBlocks.print_annotation(io::IO, bag::Bag)
    printstyled(io, "[⊞] "; bold=true, color=:blue)
end
