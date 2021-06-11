function [decoderOut, iterationOut, Lq_out, parityCheck, LLR_out] = LDPC_Decoder(LLR, Hpos, parameters, Lq_in)
Hpos = sortrows(Hpos,1)-1;
parametersArray = [parameters.MaximumIterationCount strcmp(parameters.OutputValue, 'Whole Codeword') parameters.maxUnchangedIterations strcmp(parameters.DecisionMethod, 'Soft Decision') parameters.Threads parameters.N round((1-parameters.Rate)*parameters.N) parameters.normfactor];
if(strcmp(parameters.Algorithm, 'Flooding'))
    if(strcmp(parameters.Algorithm2,'Sum Product'))
        [decoderOut, iterationOut, Lq_out, parityCheck, LLR_out] = LDPC_Decoder_mex(parametersArray,Hpos(:,2),Hpos(:,1),LLR,Lq_in);
    else
        [decoderOut, iterationOut, Lq_out, parityCheck, LLR_out] = LDPC_Decoder_MinSum_mex(parametersArray,Hpos(:,2),Hpos(:,1),LLR,Lq_in);
    end
else
    if(strcmp(parameters.Algorithm2,'Sum Product'))
        [decoderOut, iterationOut, Lq_out, parityCheck, LLR_out] = LDPC_Decoder_Layered_mex(parametersArray,Hpos(:,2),Hpos(:,1),LLR,Lq_in);
    else
        [decoderOut, iterationOut, Lq_out, parityCheck, LLR_out] = LDPC_Decoder_Layered_MinSum_mex(parametersArray,Hpos(:,2),Hpos(:,1),LLR,Lq_in);
    end
end