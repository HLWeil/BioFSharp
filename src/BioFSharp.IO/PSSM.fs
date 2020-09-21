namespace BioFSharp.IO


open FSharpAux
open BioFSharp

open System.IO
open BioFSharp
open PositionSpecificScoringMatrices
open PositionMatrix

/// Functions for reading and writing position specific scoring matrices
module PositionSpecificScoringMatrices = 
    
    module Write = 

        /// WRITE FUNCTIONS
        let writeTransfac (sw:System.IO.StreamWriter) alphabet (nameF : 'tag -> string) (valueF : 'value -> string) (pssm:TaggedPSSM<'tag,'a,'value>) =

            let name = nameF pssm.Tag

            if Array.isEmpty alphabet then failwithf "Alphabet of matrix %s is empty" name

            sw.WriteLine("XX")
            sw.WriteLine(sprintf "DE %s" name)

            alphabet
            |> Array.map (BioItem.symbol >> string)
            |> Array.append [|"P0"|]
            |> Array.reduce (fun a b -> a + "\t" + b)
            |> sw.WriteLine

            Array.init pssm.Matrix.Length (fun i -> 
                pssm.Matrix.GetAlphabetScoresOfPosition(i)
                |> Array.map (valueF)
                |> Array.append [|sprintf "%i" (i+1)|]
                |> Array.reduce (fun a b -> a + "\t" + b)
                |> sw.WriteLine
            )
            |> ignore
            sw.WriteLine("XX")
            sw.WriteLine("//")
            sw.Flush()

        /// WRITE FUNCTIONS
        let writeJaspar (sw:System.IO.StreamWriter) alphabet (nameF : 'tag -> string) (valueF : 'value -> string) (pssm:TaggedPSSM<'tag,'a,'value>) =

            let name = nameF pssm.Tag


            if Array.isEmpty alphabet then failwithf "Alphabet of matrix %s is empty" name

            sw.WriteLine(sprintf ">%s" name)
        
            alphabet
            |> Array.iter (fun a ->
                let values = pssm.Matrix.GetPositionScoresOfItem(a) |> Array.map valueF
                let symbol = BioItem.symbol a |> string

                [|"]"|]
                |> Array.append values
                |> Array.append [|symbol;"["|]
                |> Array.reduce (fun a b -> a + "\t" + b)
                |> sw.WriteLine
            )  
            sw.Flush()

        /// Write a single pssm to a specified path, specifying the format with the writeFunc. The alphabet is taken from the matrix.
        ///
        /// nameF is used to transform the tag of the matrix into a name and valueF is used to transform the position scores to a string (e.g. for a position weight matrix one could use 'sprintf "%.3f"'
        let writeSinglePSSM writeFunc (nameF : 'tag -> string) (valueF : 'value -> string) (path:string) (pssm:TaggedPSSM<'tag,'a,'value>) =
            
            let alphabet = pssm.Matrix.Alphabet

            use sw = new System.IO.StreamWriter(path)

            writeFunc sw alphabet nameF valueF pssm

        /// Write a single pssm to a specified path, specifying the format with the writeFunc. The alphabet is specified by the user 
        ///
        /// nameF is used to transform the tag of the matrix into a name and valueF is used to transform the position scores to a string (e.g. for a position weight matrix one could use 'sprintf "%.3f"'
        let writeSinglePSSMWithAlphabet writeFunc alphabet (nameF : 'tag -> string) (valueF : 'value -> string) (path:string) (pssm:TaggedPSSM<'tag,'a,'value>) =
            
            if alphabet |> Array.isEmpty then failwith "MatrixSimilarity error: Alphabets are not allowed to be empty"

            if pssm.Matrix.Alphabet <> alphabet then printfn "MatrixSimilarity Warning: Alphabet of %s and given alphabet do not match" (nameF pssm.Tag)

            use sw = new System.IO.StreamWriter(path)

            writeFunc sw alphabet nameF valueF pssm

        /// Write multiple pssms into a single file to a specified path, specifying the format with the writeFunc. The alphabet is taken from the matrices, if they match.
        ///
        /// nameF is used to transform the tag of the matrix into a name and valueF is used to transform the position scores to a string (e.g. for a position weight matrix one could use 'sprintf "%.3f"'
        let writePSSMs writeFunc (nameF : 'tag -> string) (valueF : 'value -> string) (path:string) (pssms:TaggedPSSM<'tag,'a,'value> seq) =
            
            let alphabet = 
                pssms
                |>  Seq.fold 
                    (fun alph pssm -> 
                        let alph' = pssm.Matrix.Alphabet
                        if alph <> alph' then failwithf "MatrixSimilarity error: Alphabet of Matrix %A does not match preceeding matrix alphabets of the input array" (nameF pssm.Tag)
                        alph
                    ) 
                    (pssms |> Seq.item 0).Matrix.Alphabet

            use sw = new System.IO.StreamWriter(path)

            pssms
            |> Seq.iter (writeFunc sw alphabet nameF valueF)

        /// Write multiple pssms into a single file to a specified path, specifying the format with the writeFunc. The alphabet is specified by the user.
        ///
        /// nameF is used to transform the tag of the matrix into a name and valueF is used to transform the position scores to a string (e.g. for a position weight matrix one could use 'sprintf "%.3f"'
        let writePSSMsWithAlphabet writeFunc alphabet (nameF : 'tag -> string) (valueF : 'value -> string) (path:string) (pssms:TaggedPSSM<'tag,'a,'value> seq) =
            
            if alphabet |> Array.isEmpty then failwith "MatrixSimilarity error: Alphabets are not allowed to be empty"

            use sw = new System.IO.StreamWriter(path)

            pssms
            |> Seq.iter (fun pssm ->
                if pssm.Matrix.Alphabet <> alphabet then printfn "MatrixSimilarity Warning: Alphabet of %A and given alphabet do not match" (nameF pssm.Tag)
                writeFunc sw alphabet nameF valueF pssm
            )

    module Read = 

        open FSharpAux.IO

        let private same_group l =             
            not (String.length l = 0 || l.[0] <> '>')

        let private parseJasparMatrixRow (alphabetF : char -> #IBioItem) (valueF : string -> 'value) (line:string) =
            let sl = line  |> String.split '\t' |> Array.collect (String.split ' ') |> Array.filter ((<>) "")
            let a = sl.[0] |> char |> alphabetF

            let mutable insideBrackets = false
            let vals =            
                sl
                |> Array.skip 0
                |> Array.choose (fun x -> 
                    match x with 
                    | "[" -> 
                        insideBrackets <- true
                        None
                    | "]" -> 
                        insideBrackets <- false
                        None
                    | x when insideBrackets -> 
                        valueF x 
                        |> Some
                    | _ -> None                
                )
            a, vals
            

        let private parseJasparMatrix (alphabetF : char -> #IBioItem) (valueF : string -> 'value) (lines : string seq) : (TaggedPSSM<string,'a,'value>) =
            let tm = lines |> Seq.toArray
            if Array.length tm < 2 then raise (System.Exception "Incorrect Jaspar format")
            match tm.[0], Array.skip 1 tm with
            | (t:string), m when t.StartsWith ">" ->  
                let name = t.Remove(0,1)
                let alphabet,values = 
                    m 
                    |> Array.map (parseJasparMatrixRow alphabetF valueF)
                    |> Array.unzip
                let matrix = BaseMatrix(alphabet,values.[0].Length)
                (alphabet,values) ||> Array.map2 (fun a column ->
                    column
                    |> Array.mapi (fun i v -> matrix.[a,i] <- v)
                
                ) |> ignore
                createTaggedPSSM name matrix
                                                            
            | _ -> raise (System.Exception "Incorrect Jaspar format")        
        
        let private parseCisBPMatrix (alphabetF : char -> #IBioItem) (valueF : string -> 'value) (lines : string seq) : BaseMatrix<'a,'value> =
            let alphabet = lines |> Seq.head |> String.split '\t' |> Array.skip 1 |> Array.map (char >> alphabetF)
            let values = 
                lines
                |> Seq.toArray
                |> Array.skip 1
                |> Array.map (fun l -> 
                    l
                    |> String.split '\t' |> Array.skip 1
                    |> Array.map valueF
                )
            let matrix = BaseMatrix(alphabet,values.Length)
            values 
            |> Array.mapi (fun i column ->
                (alphabet,column)
                ||> Array.map2 (fun a v -> matrix.[a,i] <- v)
            
            ) |> ignore
            matrix


        let readSingleJasparMatrix (alphabetF : char -> #IBioItem) (valueF : string -> 'value) (path:string) =
            FileIO.readFile path
            |> Seq.filter (fun (l:string) -> (not (l.StartsWith ";" || l.StartsWith "#")) && l <> "")
            |> parseJasparMatrix alphabetF valueF

        let readMultipleJasparMatrices (alphabetF : char -> #IBioItem) (valueF : string -> 'value) (path:string) =
            FileIO.readFile path
            |> Seq.filter (fun (l:string) -> (not (l.StartsWith ";" || l.StartsWith "#")) && l <> "")            
            |> Seq.groupWhen same_group 
            |> Seq.map (parseJasparMatrix alphabetF valueF)

        let readCisBPMatrix (alphabetF : char -> #IBioItem) (valueF : string -> 'value) (name:string) (path:string) =
            FileIO.readFile path
            |> Seq.filter (fun (l:string) -> (not (l.StartsWith ";" || l.StartsWith "#")) && l <> "")
            |> parseCisBPMatrix alphabetF valueF
            |> createTaggedPSSM name
        //// main
        //fileEnumerator
        //|> Seq.filter (fun (l:string) -> not (l.StartsWith ";" || l.StartsWith "#"))
        //|> Seq.groupWhen same_group 
        //|> Seq.map (fun l -> record (List.ofSeq l) converter)
        

        ///// Reads FastaItem from file. Converter determines type of sequence by converting seq<char> -> type
        //let fromFile converter (filePath) =
        //    FileIO.readFile filePath
        //    |> fromFileEnumerator converter
