" Vim syn file
" Language: Workflow Description Language
" Maintainer: Scott Frazer
" Latest Revision: 21 July 2015
"
if exists("b:current_syn")
   finish
endif

" command { ... } section
syntax region wdlCommandSection start="command\s*{" end="\v\}" contains=wdlCommand,wdlCommandParameter,wdlKeyword,wdlCommandDelimiter
syntax region wdlCommandSection2 start="command\s*<<<" end="\v>>>" contains=wdlCommand,wdlCommandParameter,wdlKeyword,wdlCommandDelimiter
syntax keyword wdlCommandKeyword command contained containedin=wdlCommandSection
syntax match wdlCommand "\zs.\{-}\ze\${" contained containedin=wdlCommandSection
syntax region wdlCommandParameter start=/\v\$\{/ end=/\v\}/ oneline contained containedin=wdlCommandSection contains=wdlType,wdlString,wdlCommandParameterName
syntax match wdlCommandParameterName /\v\zs\w+\ze([\\?\\*\\+]?\})/ contained containedin=wdlCommandParameter

" Keywords
syntax keyword wdlKeyword workflow task call nextgroup=wdlTaskName
syntax keyword wdlKeyword output scatter if then else runtime
syntax keyword wdlType Boolean Int Float String File Uri nextgroup=wdlIdentifier
syntax keyword wdlImport import

" Compound Types
syntax region wdlType start=/\(Map\|Array\)\[/ end=/\]/ contains=wdlType nextgroup=wdlIdentifier

" Identifiers
syntax match wdlIdentifier /\v\s*\w+/ contained
syntax match wdlTaskName /\v\s*\w+/ contained

" Strings
syntax region wdlString start=/"/ skip=/\\"/ end=/"/ oneline contains=wdlInterpolationWrapper
syntax region wdlInterpolationWrapper start="\v\$\{" end="\v\}" contained containedin=wdlString contains=wdlInterpolatedString
syntax match wdlInterpolatedString "\v\w+" contained containedin=wdlInterpolationWrapper

" Comments
syntax match wdlComment "\v#.*$"

highlight link wdlCommandParameter Comment
highlight link wdlKeyword Keyword
highlight link wdlCommandKeyword Keyword
highlight link wdlCommand Punctuation
highlight link wdlTaskName Identifier

highlight link wdlCommandParameterName Identifier
highlight link wdlIdentifier Identifier
highlight link wdlType Type
highlight link wdlString String
highlight link wdlImport Include
highlight link wdlInterpolationWrapper Include
highlight link wdlInterpolatedString Include
highlight link wdlComment Comment

setlocal commentstring=//\ %s
" @-@ adds the literal @ to iskeyword for @IBAction and similar
setlocal tabstop=2
setlocal softtabstop=2
setlocal shiftwidth=2

au BufRead,BufNewFile *.wdl set filetype=wdl
