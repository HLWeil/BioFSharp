namespace BioFSharp

open System
open FSharpAux
open BioFSharp

(*
TODO

-Backgroundfrequencies:
    rename \
    make functions more concise
    add table for model species

Similarity \
    1vsMany \
    ManyVsEachOther \

IO
    Read PWMS
    Write PWMS
    Hier oder in BioFSharp.IO? \

Dokumentation

*)
/// Functions for creating, comparing and using position specific scoring matrices (PWeightM, PFrequencyM, PPropabilityM).
module PositionSpecificScoringMatrices = 

    type IBioItemFrequencies =
        abstract member Print: unit -> string
    

    ///Consists of functions to work with Frequency- and Probability-Composite-Vectors.
    module BackgroundFrequencies =

        /// One dimensional array with fixed positions for each element.
        type CompositeVector<'a, 'value when 'a :> IBioItem>internal (alphabet:'a [],array:'value []) =
    
            let mutable internalAlphabet = alphabet

            let getIndex (key:'a) =
                (BioItem.symbol key |> int) - 42

            new () =
                let arr:'value [] = Array.zeroCreate 49
                new CompositeVector<_, 'value>(Array.empty,arr)

            member this.Array = 
                if (obj.ReferenceEquals(array, null)) then
                    raise (ArgumentNullException("array"))
                array
            
            member this.Item
                with get i       = array.[getIndex i]
                and  set i value = array.[getIndex i] <- value
                        
            member this.Alphabet
                with get ()                     = internalAlphabet
                and set (alphabet)  = internalAlphabet <- alphabet


        let private printHelper (toString : 'value -> string) (baseVector : CompositeVector<'a,'value>) =
            if baseVector.Alphabet |> Array.isEmpty then
                sprintf "%A" baseVector.Array
            else
                let symbols,values = 
                    baseVector.Alphabet
                    |> Array.map (fun item ->
                        let s = 
                            baseVector.Item(item)
                            |> toString    
                        if s.Length >= 5 then s,sprintf "%c\t" item.Symbol
                        else s,sprintf "%c" item.Symbol
                    )
                    |> Array.unzip 
                    |> fun (vs,ss) -> 
                        ss |> Array.reduce (fun a b -> a + "\t" + b),
                        vs |> Array.reduce (fun a b -> a + "\t" + b)
                symbols + "\n" + values

        /// One dimensional array with fixed positions for each element.
        /// Use to track frequency of elements independent of position in source.
        type FrequencyCompositeVector<'a when 'a :> IBioItem> internal (alphabet:'a [],array:int []) =

            inherit CompositeVector<'a, int>(alphabet,array)

            new() = 
                let arr = Array.zeroCreate 49
                new FrequencyCompositeVector<'a>(Array.empty,arr)

            new(alphabet) = 
                let arr = Array.zeroCreate 49
                new FrequencyCompositeVector<'a>(alphabet,arr)

            interface IBioItemFrequencies with
                member this.Print() = printHelper (sprintf "%i") this

        /// Increase counter of position by 1.
        let internal increaseInPlaceFCV (bioItem:#IBioItem) (frequencyCompositeVector:FrequencyCompositeVector<'a>) = 
            frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + 1
            frequencyCompositeVector

        /// Increase counter of position by n.
        let internal increaseInPlaceFCVBy (bioItem:'a) (frequencyCompositeVector:FrequencyCompositeVector<'a>) n = 
            frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + n
            frequencyCompositeVector

        /// Create a FrequencyCompositeVector based on BioArrays and exclude the specified segments.
        let createFCVOfSequence (resSources:BioArray.BioArray<#IBioItem>) =
            let alphabet = resSources |> Array.distinct
            let backGroundCounts = new FrequencyCompositeVector<'a>(alphabet)   
            Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts resSources

        /// Create new FrequencyCompositeVector based on the sum of the positions of an array of FrequencyVectors.
        let fuseFrequencyVectors (alphabet:#IBioItem[]) (bfVectors:FrequencyCompositeVector<'a>[]) =
            let backgroundFrequencyVector = new FrequencyCompositeVector<'a>(alphabet)
            for fcVector in bfVectors do
                for bioItem in alphabet do
                    backgroundFrequencyVector.[bioItem] <- backgroundFrequencyVector.[bioItem] + (fcVector.[bioItem])
            backgroundFrequencyVector

        /// Create FrequencyCompositeVector based on BioArray and excludes the specified segment.
        let createFCVWithout (motiveLength:int) (position:int) (resSource:BioArray.BioArray<'a>) =
            let backGroundCounts = new FrequencyCompositeVector<'a>()
            Array.append resSource.[0..(position-1)] resSource.[(position+motiveLength)..]
            |> Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts

        /// Create FrequencyCompositeVector based on BioArray.
        let internal increaseInPlaceFCVOf (resSources:BioArray.BioArray<'a>) (backGroundCounts:FrequencyCompositeVector<'a>) =   
            resSources
            |> Array.fold (fun bc bioItem -> (increaseInPlaceFCV bioItem bc)) backGroundCounts

        /// Subtracts the amount of elements in the given source from the FrequencyCompositeVector.
        let substractSegmentCountsFrom (source:BioArray.BioArray<'a>) (fcVector:FrequencyCompositeVector<'a>) =
            let bfVec = new FrequencyCompositeVector<'a>(fcVector.Alphabet,fcVector.Array)
            for bioItem in source do
                bfVec.[bioItem] <- (if fcVector.[bioItem] - 1 > 0 then fcVector.[bioItem] - 1 else 0)
            bfVec
    
        type ProbabilityCompositeVector<'a when 'a :> IBioItem> internal (alphabet:'a [],array:float []) =
        
            inherit CompositeVector<'a,float>(alphabet,array)

            new() = 
                let arr = Array.zeroCreate 49
                new ProbabilityCompositeVector<'a>(Array.empty,arr)

            new(alphabet) = 
                let arr = Array.zeroCreate 49
                new ProbabilityCompositeVector<'a>(alphabet,arr)

            interface IBioItemFrequencies with
                member this.Print() = printHelper (sprintf "%.3f") this

        /// Increase counter of position by 1.
        let internal increaseInPlacePCV (bioItem:'a) (backGroundProbabilityVector:ProbabilityCompositeVector<'a>) = 
            backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + 1.
            backGroundProbabilityVector

        /// Increase counter of position by n.
        let internal increaseInPlacePCVBy (bioItem:'a) n (backGroundProbabilityVector:ProbabilityCompositeVector<'a>) = 
            backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + n
            backGroundProbabilityVector

        /// Create a ProbabilityCompositeVector based on FrequencyCompositeVector.
        let createPCVOfFCV (alphabet:#IBioItem[]) (pseudoCount:float) (frequencyCompositeVector:FrequencyCompositeVector<'a>) =
            let backGroundProbabilityVector = 
                frequencyCompositeVector.Array
                |> Array.map (fun item -> float item)
                |> fun item -> new ProbabilityCompositeVector<'a>(alphabet,item)
            let sum = (float (Array.sum frequencyCompositeVector.Array)) + ((float alphabet.Length) * pseudoCount)
            for item in alphabet do
                backGroundProbabilityVector.[item] <- (backGroundProbabilityVector.[item] + pseudoCount)/sum
            backGroundProbabilityVector

        /// Calculate the score of the given sequence based on the probabilities of the ProbabilityCompositeVector.
        let calculateSegmentScoreBy (pcv:ProbabilityCompositeVector<'a>) (bioItems:BioArray.BioArray<#IBioItem>) =
            Array.fold (fun (value:float) (bios:#IBioItem) -> value * (pcv.[bios])) 1. bioItems

    ///Consists of functions to work with Position-Frequancy-, -Probability- and -WeightMatrices.
    module PositionMatrix =     

        ///Checks whether all elements in the list have a wider distance than width or not.
        let internal checkForDistance (width:int) (items:int list) =
            if items.Length <= 0 || items.Length = 1 then true
            else
                let rec loop n i =
                    if n = items.Length-1 then true
                    else
                        if i >= items.Length then loop (n+1) (n+2)
                        else
                            if Operators.abs(items.[n]-items.[i]) > width then
                                loop n (i+1)
                            else false
                loop 0 1

        /// Get an integer which is between 0 and the length of the sequence - segmentLength
        let internal getRandomNumberInSequence (segmentLength:int) (source:'a[]) =
            let rnd = System.Random()
            Array.init 1 (fun _ -> rnd.Next(0, source.Length-segmentLength+1))
            |> Array.head

        /// Create a specific sub sequence of the source sequence based on the given length and starting position. Do not forget, sequences start counting at 0!
        let internal getSegment (subsequenceLength:int) (source:'a[]) (startPoint:int) =
            source 
            |> Array.skip startPoint 
            |> Array.take subsequenceLength 
            |> Array.ofSeq

        ///Finds the best information content in an array of arrays of PWMSs and positions.
        let getBestInformationContent (item:((float*int)[])[]) =
            let rec loop (n:int) (bestPWMS:(float*int)[]) =
                if n = item.Length then bestPWMS
                else
                    let informationContentItem =
                        item.[n]
                        |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                    if informationContentItem > informationContentBestPWMS then
                        loop (n + 1) item.[n]
                    else
                        loop (n + 1) bestPWMS
            loop 0 [||]

        /// Matrix with fixed positions for nucleotides and amino acids.
        type BaseMatrix<[<EqualityConditionalOn; ComparisonConditionalOn >]'a, 'value when 'a :> IBioItem> internal (alphabet:'a [],matrix:'value [,]) =

            let mutable internalAlphabet = alphabet

            let getRowArray2DIndex (key:'a) =
                (BioItem.symbol key |> int) - 42

            new (rowLength:int) =
                let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
                new BaseMatrix<_, 'value>(Array.empty,arr)

            new (alphabet,rowLength:int) =
                let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
                new BaseMatrix<_, 'value>(alphabet,arr)

            member this.Matrix = 
                if (obj.ReferenceEquals(matrix, null)) then
                    raise (ArgumentNullException("array"))
                matrix
            
            member this.Item
                with get (column, row)       = matrix.[getRowArray2DIndex column, row]
                and  set (column, row) value = matrix.[getRowArray2DIndex column, row] <- value            

            member this.Length =
                matrix |> Array2D.length2

            member this.Alphabet
                with get ()                     = internalAlphabet
                and set (alphabet)  = internalAlphabet <- alphabet
              
            member this.GetPositionScoresOfItem(item) =
                Array.init this.Length (fun i ->
                    this.Item(item,i)
                )

            member this.GetAlphabetScoresOfPosition(position) =
                Array.map (fun item ->
                    this.Item(item,position)
                ) this.Alphabet

        /// Position specific scoring matrix with tag
        type TaggedPSSM<'tag,'a,'value when 'a :> IBioItem> = 
            {
                Tag     : 'tag
                Matrix  : BaseMatrix<'a,'value>
            }


        /// Auxiliary function for creating tagged PSSM
        let createTaggedPSSM tag matrix =
            {
                Tag     = tag
                Matrix  = matrix
            }

        let private printHelper (toString : 'value -> string) (baseMatrix : BaseMatrix<'a,'value>) =
            if baseMatrix.Alphabet |> Array.isEmpty then
                sprintf "%A" baseMatrix.Matrix
            else                       
                baseMatrix.Alphabet
                |> Array.map (fun item ->
                    baseMatrix.GetPositionScoresOfItem(item)
                    |> Array.fold (fun state value -> state + " " + toString value) ""
                    |> sprintf "%c %s" item.Symbol                        
                )
                |> Array.reduce (fun a b -> a + "\n" + b)


        /// Matrix representing a sequence pattern, where for each position the number of observations of a bioitem at the given position is stored
        type PositionFrequencyMatrix<'a when 'a :> IBioItem> internal (alphabet,matrix:int [,]) =

            inherit BaseMatrix<'a, int>(alphabet,matrix)

            new (rowLength:int) = 
                let arr = Array2D.zeroCreate 49 rowLength 
                new PositionFrequencyMatrix<'a>(Array.empty,arr)

            new (alphabet,rowLength:int) = 
                let arr = Array2D.zeroCreate 49 rowLength 
                new PositionFrequencyMatrix<'a>(alphabet,arr)

            interface IBioItemFrequencies with
                member this.Print() = 
                    printHelper (sprintf "%i") this

        /// Increase counter of PositionFrequencyMatrix at fixed position by 1.
        let internal increaseInPlacePFM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionFrequencyMatrix:PositionFrequencyMatrix<'a>) = 
            positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + 1
            positionFrequencyMatrix

        /// Increase counter of PositionFrequencyMatrix at fixed position by n.
        let internal increaseInPlacePFMBy (pos:int) (bioItem:'a when 'a :> IBioItem) n (positionFrequencyMatrix:PositionFrequencyMatrix<'a>) = 
            positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + n
            positionFrequencyMatrix

        /// Create PositionFrequencyMatrix based on BioArray.
        let createPFMOfSingleSequence (source:BioArray.BioArray<#IBioItem>) =
            let positionFrequencyMatrix = new PositionFrequencyMatrix<#IBioItem>(source.Length)
            source
            |> Array.fold (fun (row, cm) column -> row + 1, increaseInPlacePFM row column cm) (0, positionFrequencyMatrix) |> ignore
            positionFrequencyMatrix

        /// Create new PositionFrequencyMatrix based on the sum of the positions of an array of Countmatrices.
        let fusePositionFrequencyMatrices (motiveLength:int) (countMatrices:PositionFrequencyMatrix<'a>[]) =
            let positionFrequencyMatrix = new PositionFrequencyMatrix<'a>(motiveLength)
            if countMatrices.Length > 0 then
                for cMatrix in countMatrices do
                    for column = 0 to (Array2D.length1 cMatrix.Matrix)-1 do
                            for row = 0 to (Array2D.length2 cMatrix.Matrix)-1 do
                                positionFrequencyMatrix.Matrix.[column, row] <- positionFrequencyMatrix.Matrix.[column, row] + (cMatrix.Matrix.[column, row])
                positionFrequencyMatrix
            else positionFrequencyMatrix 
        
        //let createPFMOfMultipleSequences

        /// Matrix representing a sequence pattern, where for each position the probability to observe a specific bioitem is stored
        type PositionProbabilityMatrix<'a when 'a :> IBioItem> internal (alphabet,matrix:float [,]) =

            inherit BaseMatrix<'a, float>(alphabet,matrix)

            new(rowLength:int) = 
                let arr = Array2D.zeroCreate 49 rowLength 
                new PositionProbabilityMatrix<'a>(Array.empty,arr)

            new(alphabet,rowLength:int) = 
                let arr = Array2D.zeroCreate 49 rowLength 
                new PositionProbabilityMatrix<'a>(alphabet,arr)
            
            interface IBioItemFrequencies with
                member this.Print() = 
                    printHelper (sprintf "%.3f") this


        /// Increase counter of PositionProbabilityMatrix at fixed position by 1.
        let internal increaseInPlacePPM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionProbabilityMatrix:PositionProbabilityMatrix<'a>) = 
            positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + 1.
            positionProbabilityMatrix

        /// Increase counter of PositionProbabilityMatrix at fixed position by n.
        let internal increaseInPlacePPMBy (pos:int) (bioItem:'a when 'a :> IBioItem) n (positionProbabilityMatrix:PositionProbabilityMatrix<'a>) = 
            positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + n
            positionProbabilityMatrix

        /// Create new PositionWeightMatrix based on existing PositionFrequencyMatrix. 
        /// The frequencies are increased by the pseudoCount and for each position, the probability against all other symbols given by the alphabet is calculated.
        let createPPMOfPFM (sourceCount:int) (alphabet:#IBioItem[]) (pseudoCount:float) (positionFrequencyMatrix:PositionFrequencyMatrix<'a>) =
            let positionProbabilityMatrix = 
                positionFrequencyMatrix.Matrix |> Array2D.map (fun item -> float item)
                |> fun item -> new PositionProbabilityMatrix<'a>(alphabet,item)
            let sum = (float sourceCount) + ((float alphabet.Length) * pseudoCount)
            for item in alphabet do
                for position = 0 to (Array2D.length2 positionProbabilityMatrix.Matrix) - 1 do
                positionProbabilityMatrix.[item, position] <- (positionProbabilityMatrix.[item, position] + pseudoCount)/sum
            positionProbabilityMatrix
    
        /// Matrix representing a sequence pattern, where for each position the importance of a specific bioitem is scored against the background.
        type PositionWeightMatrix<'a when 'a :> IBioItem> internal (alphabet,matrix:float [,]) =

            inherit BaseMatrix<'a, float>(alphabet,matrix)

            new(rowLength:int) = 
                let arr = Array2D.zeroCreate 49 rowLength 
                new PositionWeightMatrix<'a>(Array.empty,arr)

            new(alphabet,rowLength:int) = 
                let arr = Array2D.zeroCreate 49 rowLength 
                new PositionWeightMatrix<'a>(alphabet,arr)

            interface IBioItemFrequencies with
                member this.Print() = 
                    printHelper (sprintf "%.3f") this

        
    
        /// Increase counter of PositionWeightMatrix at fixed position by 1.
        let internal increaseInPlacePWM (pos:int) (bioItem:'a when 'a :> IBioItem) (positionWeightMatrix:PositionWeightMatrix<'a>) = 
            positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + 1.
            positionWeightMatrix

        // Increase counter of PositionWeightMatrix at fixed position by n.
        let internal increaseInPlacePWMBy (bioItem:'a when 'a :> IBioItem) (pos:int) n (positionWeightMatrix:PositionWeightMatrix<'a>) = 
            positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + n
            positionWeightMatrix

        /// Create PositionWeightMatrix based on ProbabilityCompositeVector and PositionProbabilityMatrix.
        let createPositionWeightMatrix (alphabet:#IBioItem[]) (pcv:BackgroundFrequencies.ProbabilityCompositeVector<'a>) (ppMatrix:PositionProbabilityMatrix<'a>) =
            let pwMatrix = new PositionWeightMatrix<'a>(Array2D.length2 ppMatrix.Matrix)
            pwMatrix.Alphabet <- alphabet
            for item in alphabet do
                for position=0 to (Array2D.length2 ppMatrix.Matrix)-1 do
                    pwMatrix.[item, position] <- ppMatrix.[item, position]/pcv.[item]
            pwMatrix

        /// Calculate the score of the given sequence based on the Probabilities of the PositionWeightMatrix.
        let calculateSegmentScoreBy (pwMatrix:PositionWeightMatrix<'a>) (bioItems:BioArray.BioArray<#IBioItem>) =
            Array.fold (fun (position:int, value:float) (bios:#IBioItem) -> 
                position + 1, value * (pwMatrix.[bios, position])) (0, 1.) bioItems
            |> snd


        /// Casts a positionMatrix to a position frequency matrix. 
        let inline positionFrequencyMatrix (m : #BaseMatrix<'a,'value>) =
            let newPFM = PositionProbabilityMatrix(m.Alphabet,m.Length)
            m.Alphabet
            |> Array.iter (fun item ->
                newPFM.GetPositionScoresOfItem(item)
                |> Array.iteri (fun i v -> 
                    newPFM.Item(item,i) <- float v
                )            
            )    
            newPFM

        /// Casts a positionMatrix to a position frequency matrix by applying a given function to each value in the matrix. 
        let positionFrequencyMatrixWith (f : 'value -> int) (m : #BaseMatrix<'a,'value>) =
            let newM = m.Matrix |> Array2D.map f
            PositionFrequencyMatrix(m.Alphabet,newM)

        /// Casts a position Matrix to a position probability matrix. 
        let inline positionProbabilityMatrix (m : #BaseMatrix<'a,'value>) =
            let newPPM = PositionProbabilityMatrix(m.Alphabet,m.Length)
            m.Alphabet
            |> Array.iter (fun item ->
                newPPM.GetPositionScoresOfItem(item)
                |> Array.iteri (fun i v -> 
                    newPPM.Item(item,i) <- float v
                )            
            )            
            newPPM

        /// Casts a position Matrix to a position probability matrix by applying a given function to each value in the matrix. 
        let positionProbabilityMatrixWith (f : 'value -> float) (m : #BaseMatrix<'a,'value>) =
            let newM = m.Matrix |> Array2D.map f
            PositionProbabilityMatrix(m.Alphabet,newM)

        /// Casts a position Matrix to a position weight matrix. 
        let inline positionWeightMatrix (m : #BaseMatrix<'a,'value>) =
            let newPWM = PositionProbabilityMatrix(m.Alphabet,m.Length)
            m.Alphabet
            |> Array.iter (fun item ->
                newPWM.GetPositionScoresOfItem(item)
                |> Array.iteri (fun i v -> 
                    newPWM.Item(item,i) <- float v
                )            
            )            
            newPWM
    
        /// Casts a position Matrix to a position weight matrix by applying a given function to each value in the matrix. 
        let positionWeightMatrixWith (f : 'value -> float) (m : #BaseMatrix<'a,'value>) =
            let newM = m.Matrix |> Array2D.map f
            PositionWeightMatrix(m.Alphabet,newM)


    module SiteSampler =

        open BackgroundFrequencies

        /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
        let getBestPWMSsWithBPV (motiveLength:int) (alphabet) (source:BioArray.BioArray<'a>) (pcv:ProbabilityCompositeVector<'a>) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix<'a>) =
            let rec loop (n:int) (highValue:float) (highIndex:int) =
                if n + motiveLength > source.Length then log2(highValue), highIndex
                else
                    let tmp =
                        let segment =
                            Array.skip n source
                            |> Array.take motiveLength
                        let pwMatrix = PositionMatrix.createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                        segment
                        |> PositionMatrix.calculateSegmentScoreBy pwMatrix
                    if tmp > highValue then loop (n + 1) tmp n
                    else loop (n + 1) highValue highIndex
            loop 0 0. 0

        /// Checks whether downstream of given positions a higher InformationContent is present or not. 
        /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
        let getRightShiftedBestPWMSsWithBPV rnd (motiveLength:int) (pseudoCount:float) (alphabet:'a[]) (sources:BioArray.BioArray<'a>[]) (pcv:ProbabilityCompositeVector<'a>) (startPositions:(float*int)[]) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
                if n = sources.Length then 
                    if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                    else loop 0 acc (Array.copy acc)
                else
                    let unChosenStartPositions =
                        let tmp = Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                        Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                        |> Array.map2 (fun (source:BioArray.BioArray<'a>) (_, position) -> if position <= source.Length - motiveLength - 1 then position + 1 else position) tmp
                    let unChosenArrays =
                        Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            (PositionMatrix.getSegment motiveLength subSequence position) 
                             |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays unChosenStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                    loop 
                        (n + 1) 
                        (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                            acc.[randomSourceNumber.[n]] <- tmp
                            acc
                         else acc                    
                        )
                        bestmotive
            loop 0 (Array.copy startPositions) (Array.copy startPositions)

        /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
        /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
        let getLeftShiftedBestPWMSsWithBPV rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<'a>[]) (pcv:ProbabilityCompositeVector<'a>) (startPositions:(float*int)[]) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
                if n = sources.Length then
                    if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                    else loop 0 acc (Array.copy acc)
                else
                    let unChosenStartPositions =
                        Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                        |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                    let unChosenArrays =
                        Array.append sources.[0..randomSourceNumber.[randomSourceNumber.[n]]-1] sources.[randomSourceNumber.[randomSourceNumber.[n]]+1..]
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            (PositionMatrix.getSegment motiveLength subSequence position) 
                             |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays unChosenStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                    loop 
                        (n + 1) 
                        (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                            acc.[randomSourceNumber.[n]] <- tmp
                            acc
                         else acc                    
                        )
                        bestmotive
            loop 0 (Array.copy startPositions) (Array.copy startPositions)

        /// Checks the given Sequence for the existence of a conserved motive, by scoring each segment based on the given start positions.
        /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence.
        let findBestmotiveWithStartPosition rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<'a>[]) (pcv:ProbabilityCompositeVector<'a>) (startPositions:(float*int)[]) =        
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) acc bestmotive =
                if n = sources.Length then 
                    if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                    else loop 0 acc (Array.copy acc)
                else
                    let unChosenStartPositions =
                        Array.append acc.[0..randomSourceNumber.[n]-1] acc.[randomSourceNumber.[n]+1..]
                        |> Array.map (fun (_, position) -> position)
                    let unChosenArrays =
                        Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            (PositionMatrix.getSegment motiveLength subSequence position) 
                             |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays unChosenStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    let tmp = getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix
                    loop 
                        (n + 1) 
                        (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                            acc.[randomSourceNumber.[n]] <- tmp
                            acc
                         else acc                    
                        )
                        bestmotive
            loop 0 (Array.copy startPositions) (Array.copy startPositions)

        /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
        /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
        let getPWMOfRandomStartsWithBPV rnd (motiveLength:int) (pseudoCount:float) (alphabet:'a[]) (sources:BioArray.BioArray<'a>[]) (pcv:ProbabilityCompositeVector<'a>) =    
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) acc =
                if n = sources.Length then (List.rev acc) |> List.toArray
                else
                    let unChosenArrays =
                        Array.append (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                    let randomStartPositions = 
                        unChosenArrays
                        |> Array.map (fun unChosen ->
                            PositionMatrix.getRandomNumberInSequence motiveLength unChosen)
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            (PositionMatrix.getSegment motiveLength subSequence position) 
                             |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays randomStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    loop (n + 1) (getBestPWMSsWithBPV motiveLength alphabet sources.[randomSourceNumber.[n]] pcv positionProbabilityMatrix::acc)
            loop 0 []

        /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
        /// Each iteration it is checked if better InformationContent is found or not.
        let getmotifsWithBestInformationContentWithBPV rnd (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:'a[]) (sources:BioArray.BioArray<'a>[]) (pcv:ProbabilityCompositeVector<'a>) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
                if n > numberOfRepetitions then
                    bestPWMS
                else
                    if acc = bestPWMS then 
                        bestPWMS
                    else 
                        let informationContentAcc =
                            acc
                            |> Array.map (fun (pwms, _) -> pwms)
                            |> Array.sum
                        let informationContentBestPWMS =
                            bestPWMS
                            |> Array.map (fun (pwms, _) -> pwms)
                            |> Array.sum
                        if informationContentAcc > informationContentBestPWMS then
                            loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                        else
                            let pwms =
                                getPWMOfRandomStartsWithBPV rnd motiveLength pseudoCount alphabet sources pcv
                                |> findBestmotiveWithStartPosition rnd motiveLength pseudoCount alphabet sources pcv
                                |> getLeftShiftedBestPWMSsWithBPV rnd motiveLength pseudoCount alphabet sources pcv
                                |> getRightShiftedBestPWMSsWithBPV rnd motiveLength pseudoCount alphabet sources pcv
                            loop (n + 1) pwms bestPWMS
            loop 0 [||] [|0., 0|]

        /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
        let getBestPWMSs (motiveLength:int) (alphabet:'a[]) (pseudoCount:float) (source:BioArray.BioArray<'a>) (fcVector:FrequencyCompositeVector<'a>) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix<'a>) =
            let rec loop (n:int) (highValue:float) (highIndex:int) =
                if n + motiveLength > source.Length then log2(highValue), highIndex
                else
                    let tmp =
                        let segment =
                            Array.skip n source
                            |> Array.take motiveLength
                        let pcv =
                            increaseInPlaceFCVOf source fcVector
                            |> substractSegmentCountsFrom segment 
                            |> createPCVOfFCV alphabet pseudoCount
                        let pwMatrix = PositionMatrix.createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                        segment
                        |> PositionMatrix.calculateSegmentScoreBy pwMatrix
                    if tmp > highValue then loop (n + 1) tmp n
                    else loop (n + 1) highValue highIndex
            loop 0 0. 0

        /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
        let getBestPWMSsWithPWM motiveLength (source:BioArray.BioArray<'a>) (pwMatrix:PositionMatrix.PositionWeightMatrix<'a>) =
            let rec loop (n:int) (highValue:float) (highIndex:int) =
                if n + motiveLength > source.Length then log2(highValue), highIndex
                else
                    let tmp =
                        let segment =
                            Array.skip n source
                            |> Array.take motiveLength
                        segment
                        |> PositionMatrix.calculateSegmentScoreBy pwMatrix
                    if tmp > highValue then loop (n + 1) tmp n
                    else loop (n + 1) highValue highIndex
            loop 0 0. 0

        /// Checks whether downstream of given positions a higher InformationContent is present or not. 
        /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
        let getRightShiftedBestPWMSs rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
                if n = sources.Length then 
                    if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                    else loop 0 acc (Array.copy acc)
                else
                    let unChosenStartPositions =
                        let tmp = Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                        Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                        |> Array.map2 (fun (source:BioArray.BioArray<#IBioItem>) (_, position) -> if position <= source.Length - motiveLength - 1 then position + 1 else position) tmp
                    let unChosenArrays =
                        Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                    let frequencyCompositeVector =
                        Array.map2 (fun unchosenArray position ->
                            createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                        |> fuseFrequencyVectors alphabet
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays unChosenStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                    loop 
                        (n + 1) 
                        (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                            acc.[randomSourceNumber.[n]] <- tmp
                            acc
                         else acc                    
                        )
                        bestmotive
            loop 0 (Array.copy startPositions) (Array.copy startPositions)

        /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
        /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
        let getLeftShiftedBestPWMSs rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) (acc:(float*int)[]) (bestmotive:(float*int)[]) =
                if n = sources.Length then
                    if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                    else loop 0 acc (Array.copy acc)
                else
                    let unChosenStartPositions =
                        Array.append bestmotive.[0..randomSourceNumber.[n]-1] bestmotive.[randomSourceNumber.[n]+1..]
                        |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                    let unChosenArrays =
                        Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                    let frequencyCompositeVector =
                        Array.map2 (fun unchosenArray position ->
                            createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                        |> fuseFrequencyVectors alphabet
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays unChosenStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                    loop 
                        (n + 1) 
                        (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                            acc.[randomSourceNumber.[n]] <- tmp
                            acc
                         else acc                    
                        )
                        bestmotive
            loop 0 (Array.copy startPositions) (Array.copy startPositions)

        /// Checks the given sequence for the existence of a conserved motive, by scoring each segment based on the given start positions.
        /// The new PositionWeightMatrix is calculated and updated at each iteration if segments with higher scores are found until convergence.
        let getBestPWMSsWithStartPositions rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =        
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) acc bestmotive =
                if n = sources.Length then 
                    if (acc |> Array.map (fun item -> snd item)) = (bestmotive |> Array.map (fun item -> snd item)) then acc
                    else loop 0 acc (Array.copy acc)
                else
                    let unChosenStartPositions =
                        Array.append acc.[0..randomSourceNumber.[n]-1] acc.[randomSourceNumber.[n]+1..]
                        |> Array.map (fun (_, position) -> position)
                    let unChosenArrays =
                        Array.append sources.[0..randomSourceNumber.[n]-1] sources.[randomSourceNumber.[n]+1..]
                    let frequencyCompositeVector =
                        Array.map2 (fun unchosenArray position ->
                            createFCVWithout motiveLength position unchosenArray) unChosenArrays unChosenStartPositions
                        |> fuseFrequencyVectors alphabet
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays unChosenStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    let tmp = getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                    loop 
                        (n + 1) 
                        (if fst tmp > fst acc.[randomSourceNumber.[n]] then 
                            acc.[randomSourceNumber.[n]] <- tmp
                            acc
                         else acc                    
                        )
                        bestmotive
            loop 0 (Array.copy startPositions) (Array.copy startPositions)

        /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
        /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
        let getPWMOfRandomStarts rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) (acc:((float*int)[])) =
                if n = sources.Length then acc
                else
                    let unChosenArrays =
                        Array.append (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                    let randomStartPositions = 
                        unChosenArrays
                        |> Array.map (fun unChosen ->
                            PositionMatrix.getRandomNumberInSequence motiveLength unChosen)
                    let frequencyCompositeVector =
                        Array.map2 (fun unchosenArray position ->
                            createFCVWithout motiveLength position unchosenArray) unChosenArrays randomStartPositions
                        |> fuseFrequencyVectors alphabet
                    let positionProbabilityMatrix =
                        Array.map2 (fun subSequence position -> 
                            PositionMatrix.getSegment motiveLength subSequence position
                            |> PositionMatrix.createPFMOfSingleSequence) unChosenArrays randomStartPositions
                        |> PositionMatrix.fusePositionFrequencyMatrices motiveLength
                        |> PositionMatrix.createPPMOfPFM (sources.Length - 1) alphabet pseudoCount
                    loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                                  acc)
            loop 0 (Array.zeroCreate sources.Length)

        /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
        /// Each iteration it is checked if better InformationContent is found or not.
        let getmotifsWithBestInformationContent rnd (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
            let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
                if n > numberOfRepetitions then
                    bestPWMS
                else
                    if acc = bestPWMS then 
                        bestPWMS
                    else 
                        let informationContentAcc =
                            acc
                            |> Array.map (fun (pwms, _) -> pwms)
                            |> Array.sum
                        let informationContentBestPWMS =
                            bestPWMS
                            |> Array.map (fun (pwms, _) -> pwms)
                            |> Array.sum
                        if informationContentAcc > informationContentBestPWMS then
                            loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                        else
                            let pwms =
                                getPWMOfRandomStarts rnd motiveLength pseudoCount alphabet sources
                                |> getBestPWMSsWithStartPositions rnd motiveLength pseudoCount alphabet sources
                                |> getLeftShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources
                                |> getRightShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources
                            loop (n + 1) pwms bestPWMS
            loop 0 [||] [|0., 0|]

        /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
        /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
        let getmotifsWithBestPWMSOfPPM rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<'a>[]) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix<'a>) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) acc:((float*int)[]) =
                if n = sources.Length then acc
                else
                    let unChosenArrays =
                        Array.append (sources.[0..randomSourceNumber.[n]-1]) (sources.[randomSourceNumber.[n]+1..])
                    let randomStartPositions = 
                        unChosenArrays
                        |> Array.map (fun unChosen ->
                            PositionMatrix.getRandomNumberInSequence motiveLength unChosen)
                    let frequencyCompositeVector =
                        Array.map2 (fun unchosenArray position ->
                            createFCVWithout motiveLength position unchosenArray) 
                                unChosenArrays randomStartPositions
                        |> fuseFrequencyVectors alphabet
                    loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSs motiveLength alphabet pseudoCount sources.[randomSourceNumber.[n]] frequencyCompositeVector positionProbabilityMatrix
                                  acc)
            loop 0 (Array.zeroCreate sources.Length)

        /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
        /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
        let getmotifsWithBestPWMSOfPWM rnd (motiveLength:int) (sources:BioArray.BioArray<'a>[]) (pwm:PositionMatrix.PositionWeightMatrix<'a>) =
            let randomSourceNumber = Array.shuffleFisherYates rnd ([|0..sources.Length-1|])
            let rec loop (n:int) acc:((float*int)[]) =
                if n = sources.Length then acc
                else
                    loop (n + 1) (acc.[randomSourceNumber.[n]] <- getBestPWMSsWithPWM motiveLength sources.[randomSourceNumber.[n]] pwm
                                  acc)
            loop 0 (Array.zeroCreate sources.Length)

        /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
        /// Each iteration it is checked if better InformationContent is found or not.
        let getBestInformationContentOfPPM rnd (numberOfRepetitions:int) (motiveLength:int) (pseudoCount:float) (alphabet:'a[]) (sources:BioArray.BioArray<'a>[]) (positionProbabilityMatrix:PositionMatrix.PositionProbabilityMatrix<'a>) =
            let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
                if n > numberOfRepetitions then
                    bestPWMS
                else
                    if acc = bestPWMS then 
                        bestPWMS
                    else 
                        let informationContentAcc =
                            acc
                            |> Array.map (fun (pwms, _) -> pwms)
                            |> Array.sum
                        let informationContentBestPWMS =
                            bestPWMS
                            |> Array.map (fun (pwms, _) -> pwms)
                            |> Array.sum
                        if informationContentAcc > informationContentBestPWMS then
                            loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                        else
                            let pwms =
                                getmotifsWithBestPWMSOfPPM rnd motiveLength pseudoCount alphabet sources positionProbabilityMatrix
                                |> getBestPWMSsWithStartPositions rnd motiveLength pseudoCount alphabet sources
                                |> getLeftShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources
                                |> getRightShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources
                            loop (n + 1) pwms bestPWMS
            loop 0 [||] [|0., 0|]

        let doSiteSampling rnd (motiveLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
            getPWMOfRandomStarts rnd motiveLength pseudoCount alphabet sources
            |> getBestPWMSsWithStartPositions rnd motiveLength pseudoCount alphabet sources
            |> getLeftShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources
            |> getRightShiftedBestPWMSs rnd  motiveLength pseudoCount alphabet sources

        let doSiteSamplingWithPCV rnd (motiveLength:int) (pseudoCount:float) (alphabet:'a[]) (sources:BioArray.BioArray<'a>[]) (pcv:ProbabilityCompositeVector<'a>) =
            getPWMOfRandomStartsWithBPV rnd motiveLength pseudoCount alphabet sources pcv
            |> findBestmotiveWithStartPosition rnd motiveLength pseudoCount alphabet sources pcv
            |> getLeftShiftedBestPWMSsWithBPV rnd motiveLength pseudoCount alphabet sources pcv
            |> getRightShiftedBestPWMSsWithBPV rnd motiveLength pseudoCount alphabet sources pcv

        let doSiteSamplingWithPPM rnd (motiveLength:int) (pseudoCount:float) (alphabet:'a[]) (sources:BioArray.BioArray<'a>[]) (ppM:PositionMatrix.PositionProbabilityMatrix<'a>) =
            getmotifsWithBestPWMSOfPPM rnd motiveLength pseudoCount alphabet sources ppM
            |> getBestPWMSsWithStartPositions rnd motiveLength pseudoCount alphabet sources
            |> getLeftShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources
            |> getRightShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources

        let doSiteSamplingWithPWM rnd (motiveLength:int) (pseudoCount:float) (alphabet:'a[]) (sources:BioArray.BioArray<'a>[]) (pwm:PositionMatrix.PositionWeightMatrix<'a>) =
            getmotifsWithBestPWMSOfPWM rnd motiveLength sources pwm
            |> getBestPWMSsWithStartPositions rnd motiveLength pseudoCount alphabet sources
            |> getLeftShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources
            |> getRightShiftedBestPWMSs rnd motiveLength pseudoCount alphabet sources

    /// Functions for computing the similarity between position specific scoring matrices (motif similarity)
    module Similarity = 

        open PositionMatrix

        /// Extracts only the columns of the matrix specified by the alphabet 
        let private positionMatrixToDensePositionArray (alphabet:#IBioItem[]) (matrix: BaseMatrix<#IBioItem,'T>) =
            let get position bioitem = matrix.[bioitem,position]
            Array.init matrix.Length (fun position ->
                alphabet
                |> Array.map (get position)        
            )

        /// Computes similarity between two motifs by summing the position wise scores obtained by applying the comparison wise functions over the alphabets
        ///
        /// Tries all possible gapless alignments between the motifs
        let private compareDensePositionArrays (comparisonFunction : 'T [] -> 'T [] -> float) (m1 : 'T [][]) (m2 : 'T [][]) : float = 
            if m1.Length = m2.Length then
                let unnormalizedScore = 
                    Array.init m1.Length (fun i ->
                        comparisonFunction m1.[i] m2.[i]
                    )
                    |> Array.sum
                unnormalizedScore / (float m1.Length)

             else 
                let shorterMotif,longerMotif = if m1.Length < m2.Length then m1,m2 else m2,m1
                let averageMotifLength = (float m1.Length + float m2.Length) / 2.
                let lengthDifference = longerMotif.Length - shorterMotif.Length

                Array.init (lengthDifference + 1) (fun offSet -> 
                    let unnormalizedScore = 
                        Array.init shorterMotif.Length (fun i ->
                            comparisonFunction shorterMotif.[i] longerMotif.[i+offSet]
                        )
                        |> Array.sum
                    unnormalizedScore / averageMotifLength          
                )
                |> Array.max

        /// Computes similarity between two motifs by summing the position wise scores obtained by applying the comparisonFunction over given alphabet
        ///
        /// Tries all possible gapless alignments between the motifs
        let compareTwoPositionMatricesWithAlphabet (alphabet:#IBioItem[]) (comparisonFunction : 'T [] -> 'T [] -> float) (matrix1: BaseMatrix<#IBioItem,'T>) (matrix2: BaseMatrix<#IBioItem,'T>) = 

            if alphabet |> Array.isEmpty then failwith "MatrixSimilarity error: Alphabets are not allowed to be empty"
            if matrix1.Alphabet <> alphabet then printfn "MatrixSimilarity Warning: Alphabet of Matrix 1 and given alphabet do not match"
            if matrix2.Alphabet <> alphabet then printfn "MatrixSimilarity Warning: Alphabet of Matrix 2 and given alphabet do not match"

            let m1 = positionMatrixToDensePositionArray alphabet matrix1
            let m2 = positionMatrixToDensePositionArray alphabet matrix2
            compareDensePositionArrays comparisonFunction m1 m2

        /// Computes similarity between two motifs by summing the position wise scores obtained by applying the comparisonFunction over the alphabet, which is obtained from the motifs
        ///
        /// Tries all possible gapless alignments between the motifs
        let compareTwoPositionMatrices (comparisonFunction : 'T [] -> 'T [] -> float) (matrix1: BaseMatrix<#IBioItem,'T>) (matrix2: BaseMatrix<#IBioItem,'T>) =  
            
            if matrix1.Alphabet <> matrix2.Alphabet then failwith "MatrixSimilarity error: Alphabet of Matrix 1 and Matrix 2 do not match"
            if matrix1.Alphabet |> Array.isEmpty then failwith "MatrixSimilarity error: Alphabets are not allowed to be empty"

            let m1 = positionMatrixToDensePositionArray matrix1.Alphabet matrix1
            let m2 = positionMatrixToDensePositionArray matrix1.Alphabet matrix2
            compareDensePositionArrays comparisonFunction m1 m2

        /// Computes pairwise similarity between each motif in the given collection by summing the position wise scores obtained by applying the comparisonFunction over the given alphabet
        ///
        /// Tries all possible gapless alignments between the motifs
        let comparePositionMatricesWithAlphabet (alphabet:#IBioItem[]) (comparisonFunction : 'T [] -> 'T [] -> float) (matrices: TaggedPSSM<'tag,#IBioItem,'T> []) =
            if alphabet |> Array.isEmpty then failwith "MatrixSimilarity error: Alphabets are not allowed to be empty"

            let densePositionArrays = 
                matrices 
                |> Array.map (fun m -> 
                    if m.Matrix.Alphabet <> alphabet then printfn "MatrixSimilarity Warning: Alphabet of %A and given alphabet do not match" m.Tag
                    m.Tag,
                    m.Matrix |> positionMatrixToDensePositionArray alphabet
                )

            [|
                for i = 0 to densePositionArrays.Length - 1 do
                    for j = (i + 1) to densePositionArrays.Length - 1 do

                        let n1,m1 = densePositionArrays.[i]
                        let n2,m2 = densePositionArrays.[j]

                        yield (n1,n2), compareDensePositionArrays comparisonFunction m1 m2            
            |]

        /// Computes pairwise similarity between each motif in the given collection by summing the position wise scores obtained by applying the comparisonFunction over the alphabet, which is obtained from the motifs
        ///
        /// Tries all possible gapless alignments between the motifs   
        let comparePositionMatrices (comparisonFunction : 'T [] -> 'T [] -> float) (matrices: TaggedPSSM<'tag,#IBioItem,'T> []) =
            let alphabet = 
                matrices
                |> Array.fold 
                    (fun alph matrix -> 
                        let alph' = matrix.Matrix.Alphabet
                        if alph <> alph' then failwithf "MatrixSimilarity error: Alphabet of Matrix %A does not match preceeding matrix alphabets of the input array" matrix.Tag
                        alph
                    ) 
                    matrices.[0].Matrix.Alphabet
            
            comparePositionMatricesWithAlphabet alphabet comparisonFunction matrices