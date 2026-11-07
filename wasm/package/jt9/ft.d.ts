export interface DecodeOptions {
  /** 1-3, as jt9's -d (default 2) */
  depth?: number;
  /** Search range, Hz (default 100-3200) */
  freqLow?: number;
  freqHigh?: number;
  /** Rx and Tx audio frequencies, Hz (default 1500) */
  rxFreq?: number;
  txFreq?: number;
  /** a-priori search width, Hz (default 50) */
  apWidth?: number;
  /** Blank by default: only the CQ a-priori passes run */
  myCall?: string;
  hisCall?: string;
  hisGrid?: string;
  qsoProgress?: number;
  /** nexp_decode (its low 3 bits: the contest type) */
  contest?: number;
  /** a-priori decoding (default true) */
  apDecode?: boolean;
  apCqOnly?: boolean;
  /** Quicker subtraction: the same decodes, not WSJT-X's to the bit */
  fastSubtract?: boolean;
}

/** MSK144's, as jt9 -k reads them */
export interface MskDecodeOptions {
  /** The period, s: 15 or 30 (default 15), as jt9's -p */
  period?: number;
  /** Audio frequency, Hz (default 1500), as -f */
  rxFreq?: number;
  /** Hz either side of rxFreq (default 200, at most 1000), as -F */
  tolerance?: number;
  /** 1 short pings only, 2 or 3 also the averages over frames (default 3), as -d */
  depth?: number;
  /** Blank by default */
  myCall?: string;
  hisCall?: string;
  /** MSK40 short messages to myCall and hisCall, "<MYCALL HISCALL> RPT" (default false) */
  shortMessages?: boolean;
  /** SWL mode: MSK40 messages between any calls heard, "<HHHH> RPT" for an unknown hash (default false) */
  swl?: boolean;
}

export interface MskDecode {
  snr: number;
  /** Seconds into the period */
  dt: number;
  /** Audio frequency, Hz */
  freq: number;
  message: string;
  /** The line jt9 prints */
  line: string;
}

export interface MskDecodeResult {
  decodes: MskDecode[];
  /** The number of decodes */
  count: number;
  /** jt9's stdout for the period */
  text: string;
}

export interface Decode {
  snr: number;
  /** Seconds */
  dt: number;
  /** Audio frequency, Hz */
  freq: number;
  message: string;
  /** jt9's '?' */
  lowConfidence: boolean;
  /** a-priori type, 0 for none */
  ap: number;
  quality: number;
  sync: number;
  /** The line jt9 prints */
  line: string;
}

export interface DecodeResult {
  decodes: Decode[];
  /** <DecodeFinished>'s count */
  count: number;
  /** jt9's stdout for the period */
  text: string;
}

export interface Decoder {
  readonly mode: "ft8" | "ft4" | "ft2";
  /**
   * One period's 12 kHz audio from its start (FT8 15 s, FT4 7.5 s, FT2
   * 3.75 s; shorter is zero-filled). Float32Array is -1..1, Int16Array as is.
   * utc: HHMMSS.
   */
  decode(samples: Float32Array | Int16Array, utc: number, options?: DecodeOptions): DecodeResult;
  /** Frees the decoder's memory */
  delete(): void;
}

/** MSK144: each period decoded afresh, as jt9 -k decodes each file */
export interface MskDecoder {
  readonly mode: "msk144";
  /**
   * One period's 12 kHz audio from its start (options.period long; shorter
   * is zero-filled, longer is cut as jt9 cuts it). Float32Array is -1..1,
   * Int16Array as is. utc: HHMMSS.
   */
  decode(samples: Float32Array | Int16Array, utc: number, options?: MskDecodeOptions): MskDecodeResult;
  /** Frees the decoder's memory */
  delete(): void;
}

export interface FtModule {
  Decoder: {
    new (mode: "msk144"): MskDecoder;
    new (mode: "ft8" | "ft4" | "ft2"): Decoder;
  };
}

export default function createFt(moduleOverrides?: object): Promise<FtModule>;
