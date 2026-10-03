--[[
Pandoc filter that turns the conventions of manuscript.md into an Elsevier elsarticle manuscript.

The Markdown is written to read well on its own, so this filter, not the Markdown, carries the LaTeX conventions:

- the "Abstract" section becomes the abstract, and the "References" section is dropped for the BibTeX bibliography
- the italic draft note and horizontal rules are dropped
- section numbers typed in headings are dropped, since LaTeX numbers sections; unnumbered back matter stays unnumbered
- an image followed by an italic paragraph "Figure N. ..." becomes figure N with that caption, and an italic paragraph
  "Table N. ..." followed by a table becomes that table's caption
- symbols typed as text (c₁, S₀, x_N, 𝐈, ≤, →) and expressions of variable letters (k, x, Δt, 2kx/Δt) become inline
  math, so they are set like the equations
- a paragraph that continues a sentence after a displayed equation is not indented
- inline code breaks across lines, so long file paths stay inside the margins
- with the metadata flat-figures, images are referenced by file name alone, as vector PDF when one sits beside the PNG,
  for a submission folder with no subfolders, and their source paths are listed in the file named by figure-list

Run with --shift-heading-level-by=-1, which makes the level-one heading the title. Pandoc shifts the headings after
filters run, so this filter sees the sections as level two.
]]

local UNNUMBERED = { ['Code and data availability'] = true, ['Acknowledgments'] = true, ['Data availability'] = true }
local SUBSCRIPT_DIGITS = { ['₀'] = '0', ['₁'] = '1', ['₂'] = '2', ['₃'] = '3', ['₄'] = '4', ['₅'] = '5', ['₆'] = '6',
  ['₇'] = '7', ['₈'] = '8', ['₉'] = '9' }
local MAX_BLOCKS = 100000
-- symbols the text font has no glyph for, set as math wherever they are typed in text
local MATH_SYMBOLS = { ['≤'] = '\\le', ['≥'] = '\\ge', ['→'] = '\\to', ['∈'] = '\\in', ['∞'] = '\\infty',
  ['⌈'] = '\\lceil', ['⌉'] = '\\rceil' }

local flat_figures = false
local figure_list = nil
local figure_sources = {}

local function file_exists(path)
  local handle = io.open(path, 'rb')
  if handle then
    handle:close()
    return true
  end
  return false
end

local function starts_with_label(inlines, word)
  -- whether inlines begin with "<word> <number>." as in "Figure 3." or "Table 2."
  return #inlines >= 3 and inlines[1].t == 'Str' and inlines[1].text == word and inlines[2].t == 'Space'
    and inlines[3].t == 'Str' and inlines[3].text:match('^%d+%.$') ~= nil
end

local function strip_label(inlines)
  -- the inlines after "<word> <number>. ", and the number
  local number = inlines[3].text:match('^(%d+)%.$')
  local rest = pandoc.Inlines({})
  local start = (#inlines >= 4 and inlines[4].t == 'Space') and 5 or 4
  for i = start, #inlines do
    rest:insert(inlines[i])
  end
  return rest, number
end

local function italic_label(block, word)
  -- the inlines of a paragraph that is one emphasized "<word> N. ..." caption, or nil
  if block == nil or block.t ~= 'Para' or #block.content ~= 1 or block.content[1].t ~= 'Emph' then
    return nil
  end
  local inner = block.content[1].content
  if starts_with_label(inner, word) then
    return inner
  end
  return nil
end

local function flatten_image(image)
  -- reference an image by file name, preferring the vector PDF beside a PNG
  local source = image.src
  local pdf = source:gsub('%.png$', '.pdf')
  if pdf ~= source and file_exists(pdf) then
    source = pdf
  end
  table.insert(figure_sources, source)
  image.src = pandoc.path.filename(source)
  return image
end

local function caption_figure(figure, inlines)
  local caption, number = strip_label(inlines, 'Figure')
  figure.caption = { long = { pandoc.Plain(caption) } }
  figure.identifier = 'fig:' .. number
  figure.content = figure.content:walk({
    Image = function(image)
      image.attributes.width = '100%'
      if flat_figures then
        return flatten_image(image)
      end
      return image
    end,
  })
  return figure
end

local function size_columns(tbl)
  -- give each column a share of the line width that grows with its longest cell, within limits that keep every column
  -- readable; pandoc otherwise gives the columns of a wide pipe table equal shares
  local longest = {}
  for i = 1, #tbl.colspecs do
    longest[i] = 9 -- room for a header word such as Treatment
  end
  local function measure(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if longest[i] then
          longest[i] = math.max(longest[i], math.min(utf8.len(pandoc.utils.stringify(cell.contents)) or 4, 36))
        end
      end
    end
  end
  measure(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do
    measure(body.body)
  end
  local total = 0
  for i = 1, #longest do
    total = total + longest[i]
  end
  for i, spec in ipairs(tbl.colspecs) do
    tbl.colspecs[i] = { spec[1], longest[i] / total }
  end
  return tbl
end

local function caption_table(tbl, inlines)
  local caption, number = strip_label(inlines, 'Table')
  tbl.caption = { long = { pandoc.Plain(caption) } }
  tbl.identifier = 'tab:' .. number
  return size_columns(tbl)
end

local function heading_text(block)
  return pandoc.utils.stringify(block.content)
end

local function extract_sections(blocks, meta)
  -- drop the draft note and rules, lift the abstract into the metadata, and drop the hand-written references
  local kept = pandoc.Blocks({})
  local mode, mode_level = 'body', 0
  local abstract = pandoc.Blocks({})
  for i = 1, math.min(#blocks, MAX_BLOCKS) do
    local block = blocks[i]
    if block.t == 'Header' then
      local name = heading_text(block)
      if name == 'Abstract' or name == 'References' then
        mode, mode_level = name:lower(), block.level
      elseif block.level <= mode_level then
        mode, mode_level = 'body', 0
      end
    end
    local is_rule = block.t == 'HorizontalRule'
    local is_note = block.t == 'Para' and #block.content == 1 and block.content[1].t == 'Emph'
      and pandoc.utils.stringify(block.content[1]):match('^Draft manuscript') ~= nil
    if mode == 'abstract' and block.t ~= 'Header' and not is_rule and block.t ~= 'RawBlock' then
      abstract:insert(block)
    elseif mode == 'body' and not is_rule and not is_note then
      kept:insert(block)
    end
  end
  if #abstract > 0 then
    meta.abstract = pandoc.MetaBlocks(abstract)
  end
  return kept
end

local function is_display_math(block)
  return block ~= nil and block.t == 'Para' and #block.content == 1 and block.content[1].t == 'Math'
    and block.content[1].mathtype == 'DisplayMath'
end

local function continue_after_equations(blocks)
  -- a paragraph that continues a sentence after an equation, starting in lowercase, is not indented
  for i = 2, math.min(#blocks, MAX_BLOCKS) do
    local block = blocks[i]
    local first = block.t == 'Para' and block.content[1] or nil
    if is_display_math(blocks[i - 1]) and first and first.t == 'Str' and first.text:match('^%l') then
      block.content:insert(1, pandoc.RawInline('latex', '\\noindent '))
    end
  end
  return blocks
end

local function next_table(blocks, i)
  -- the index of the table after block i, passing over raw blocks such as HTML comments, or nil
  for j = i + 1, math.min(#blocks, i + 4) do
    if blocks[j].t == 'Table' then
      return j
    end
    if blocks[j].t ~= 'RawBlock' then
      return nil
    end
  end
  return nil
end

local function sized_table(tbl)
  -- a table set a size smaller than the text, or two sizes smaller when it has eight columns or more
  local size = #tbl.colspecs >= 8 and '\\footnotesize' or '\\small'
  return { pandoc.RawBlock('latex', '\\begingroup' .. size), tbl, pandoc.RawBlock('latex', '\\endgroup') }
end

local function caption_floats(blocks)
  -- attach each italic "Figure N." paragraph to the figure before it, and each "Table N." to the table after it
  local out = pandoc.Blocks({})
  local skip_to = 0
  for i = 1, math.min(#blocks, MAX_BLOCKS) do
    local block, following = blocks[i], blocks[i + 1]
    local table_at = italic_label(block, 'Table') and next_table(blocks, i) or nil
    if i <= skip_to then
      -- already placed with the caption before it
    elseif block.t == 'Figure' and italic_label(following, 'Figure') then
      out:insert(caption_figure(block, italic_label(following, 'Figure')))
      skip_to = i + 1
    elseif table_at then
      out:extend(sized_table(caption_table(blocks[table_at], italic_label(block, 'Table'))))
      skip_to = table_at
    elseif block.t == 'Table' then
      out:extend(sized_table(block))
    else
      out:insert(block)
    end
  end
  return out
end

local function math_token(text)
  -- the LaTeX of a symbol typed as text, or nil: c₁ -> c_{1}, 𝐈 -> \mathbf{I}, x_N -> x_{N}, L_min -> L_{\mathrm{min}}
  if text == '𝐈' then
    return '\\mathbf{I}'
  end
  local letter, rest = text:match('^([A-Za-z])(.+)$')
  if letter and rest then
    local digits = {}
    if not utf8.len(rest) then
      return nil
    end
    for _, code in utf8.codes(rest) do
      local digit = SUBSCRIPT_DIGITS[utf8.char(code)]
      if digit == nil then
        digits = nil
        break
      end
      table.insert(digits, digit)
    end
    if digits and #digits > 0 then
      return letter .. '_{' .. table.concat(digits) .. '}'
    end
  end
  local base, sub = text:match('^([A-Za-z])_([A-Za-z0-9]+)$')
  if base and #sub <= 3 then
    return base .. '_{' .. ((#sub > 1 and sub:match('^[a-z]+$')) and ('\\mathrm{' .. sub .. '}') or sub) .. '}'
  end
  return nil
end

-- letters that are always variables in the text, never words or units; with digits and operators they form expressions
local VARIABLES = { x = true, k = true, C = true, N = true, L = true, D = true, Q = true, q = true, v = true, n = true,
  t = true, c = true, S = true, I = true, i = true, j = true }
local GREEK = { ['Δ'] = '\\Delta ', ['ω'] = '\\omega ', ['π'] = '\\pi ' }
local SUPERSCRIPTS = { ['²'] = '^2', ['³'] = '^3' }
local OPERATORS = '0123456789()/+=.\'−-'

local function variable_expression(text)
  -- the LaTeX of a token built only of variable letters, digits, and operators, such as 2kx/Δt or (k/N)², or nil
  if not utf8.len(text) then
    return nil
  end
  local parts, has_variable = {}, false
  for _, code in utf8.codes(text) do
    local character = utf8.char(code)
    if VARIABLES[character] or GREEK[character] then
      has_variable = true
      table.insert(parts, GREEK[character] or character)
    elseif SUPERSCRIPTS[character] then
      table.insert(parts, SUPERSCRIPTS[character])
    elseif #character == 1 and OPERATORS:find(character, 1, true) then
      table.insert(parts, character)
    elseif character == '−' then
      table.insert(parts, '-')
    else
      return nil
    end
  end
  -- a letter followed by a digit is a label such as v3, not a variable (subscripts are typed as c₁)
  if text:find('%a%d') then
    return nil
  end
  -- a run of letters alone is a variable only when it is one letter or the product kx; otherwise it is a word (in, it)
  local letters_only = not text:find('[^%a]') and utf8.len(text) > 1
  if letters_only and text ~= 'kx' then
    return nil
  end
  return has_variable and table.concat(parts) or nil
end

local function typeset_token(text)
  -- one run of text as inlines: leading punctuation, a symbol set as math, and trailing punctuation, or nil
  local lead, core, trail = text:match('^([%(%[|]*)(.-)([%.,;:%)%]|]*)$')
  local tex = core and (math_token(core) or variable_expression(core))
  if tex == nil then
    return nil
  end
  local out = pandoc.Inlines({})
  if lead ~= '' then out:insert(pandoc.Str(lead)) end
  out:insert(pandoc.Math('InlineMath', tex))
  if trail ~= '' then out:insert(pandoc.Str(trail)) end
  return out
end

local function split_math_symbols(text)
  -- the runs of text between the symbols of MATH_SYMBOLS, and each symbol as inline math
  local out, buffer = pandoc.Inlines({}), {}
  if not utf8.len(text) then
    return pandoc.Inlines({ pandoc.Str(text) })
  end
  for _, code in utf8.codes(text) do
    local character = utf8.char(code)
    if MATH_SYMBOLS[character] then
      if #buffer > 0 then out:insert(pandoc.Str(table.concat(buffer))) end
      buffer = {}
      out:insert(pandoc.Math('InlineMath', MATH_SYMBOLS[character]))
    else
      table.insert(buffer, character)
    end
  end
  if #buffer > 0 then out:insert(pandoc.Str(table.concat(buffer))) end
  return out
end

local function typeset_symbols(str)
  -- set the symbols typed in a Str as math: those the text font lacks, and subscripted variables
  local pieces = split_math_symbols(str.text)
  local out, changed = pandoc.Inlines({}), #pieces ~= 1 or pieces[1].t ~= 'Str'
  for _, piece in ipairs(pieces) do
    local typeset = piece.t == 'Str' and typeset_token(piece.text) or nil
    if typeset then
      out:extend(typeset)
      changed = true
    else
      out:insert(piece)
    end
  end
  return changed and out or nil
end

local function number_headings(block)
  -- drop typed section numbers and keep back matter unnumbered
  if block.t ~= 'Header' then
    return block
  end
  if #block.content >= 2 and block.content[1].t == 'Str' and block.content[1].text:match('^%d+[%.%d]*$')
    and block.content[2].t == 'Space' then
    block.content = { table.unpack(block.content, 3) }
  end
  if UNNUMBERED[heading_text(block)] then
    block.classes:insert('unnumbered')
  end
  return block
end

local function author_names(meta)
  -- the names of the authors, for the PDF metadata, which pandoc cannot derive from author records
  local names = {}
  for _, author in ipairs(meta.author or {}) do
    table.insert(names, pandoc.utils.stringify(author.name or author))
  end
  return table.concat(names, ', ')
end

local function breakable_code(code)
  -- inline code as \\nolinkurl, which xurl lets break anywhere, so long paths do not run into the margin
  if code.text:find('[{}%%#\\]') then
    return nil
  end
  return pandoc.RawInline('latex', '\\nolinkurl{' .. code.text .. '}')
end

function Pandoc(doc)
  doc.meta['author-meta'] = author_names(doc.meta)
  flat_figures = doc.meta['flat-figures'] == true
  figure_list = doc.meta['figure-list'] and pandoc.utils.stringify(doc.meta['figure-list']) or nil
  local blocks = continue_after_equations(caption_floats(extract_sections(doc.blocks, doc.meta)))
  blocks = blocks:walk({ Header = number_headings, Str = typeset_symbols, Code = breakable_code })
  if doc.meta.abstract then
    doc.meta.abstract = pandoc.MetaBlocks(pandoc.Blocks(doc.meta.abstract):walk({ Str = typeset_symbols }))
  end
  if figure_list then
    local handle = assert(io.open(figure_list, 'w'))
    handle:write(table.concat(figure_sources, '\n'), '\n')
    handle:close()
  end
  doc.blocks = blocks
  return doc
end
