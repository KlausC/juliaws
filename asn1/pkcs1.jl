"ASN1 defintions form PKCS-1"

module PKCS1

export REG

using Register
using Syntax

import Base.(>=)

const REG = Registry()

(>=)(s::Symbol, t::RegType) = begin assign(s, t, REG); t end

OID = OBJECT_IDENTIFIER

Version = ENUMERATED( :without : 0, :with :1 )


:RSAPublicKey >= SEQUENCE(
					:modulus => [1] : INTEGER,
					:publicExponenti => [3] : INTEGER
					)

:RSAPrivateKey >= SEQUENCE(
					:version => Version,
					EXT,
					:modulus => INTEGER,
					:publicExponent => INTEGER,
					:privateExponent => INTEGER,
					:prime1 => INTEGER,
					:prime2 => INTEGER,
					:exponent1 => INTEGER,
					:exponent2 => INTEGER,
					:coefficient => INTEGER > OPTIONAL,
					EXT,
					:otherPrimeInfos => [1] : !:OtherPrimeInfos > OPTIONAL
					)

:OtherPrimeInfos >= SEQUENCE_OF(!:OtherPrimeInfo)

:OtherPrimeInfo >= SEQUENCE(
					:prime => INTEGER,
					:exponent => INTEGER,
					:coefficient => INTEGER,
					:test => !:ChoiceExample,
					)

:ChoiceExample >= CHOICE(
					 :first => [0] : !:RSAPublicKey,
					 :second => [1] : !:RSAPrivateKey
				 )

:SetExample >= SET(
			 :first => [0] : !:RSAPublicKey,
			 :second => [1] : !:RSAPrivateKey
		 )

:SetOf >= SET_OF([42] : !:RSAPublicKey)

:Loop >= SET(
			 :data => INTEGER,
			 :next => CHOICE( :a => !:Loop, :b => NULL)
		 )


end # module
