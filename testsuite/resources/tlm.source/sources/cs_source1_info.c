"{\"format\":\"Transformational debugger info\",\"version\":1,\n\"info\":{\"name\":\"source1\",\"description\":\"\"},\n\"variables\":{\n\"$cse1\":{\"comment\":\"\",\"kind\":\"variable\",\"type\":\"Real\",\"unit\":\"\",\"displayUnit\":\"\",\"source\":{\"info\":{\"file\":\"\",\"lineStart\":0,\"lineEnd\":0,\"colStart\":0,\"colEnd\":0}}},\n\"A\":{\"comment\":\"\",\"kind\":\"parameter\",\"type\":\"Real\",\"unit\":\"\",\"displayUnit\":\"\",\"source\":{\"info\":{\"file\":\"<interactive>\",\"lineStart\":3,\"lineEnd\":3,\"colStart\":3,\"colEnd\":25},\"within\":[\"Real\"]}}\n},\n\"equations\":[{\"eqIndex\":0,\"tag\":\"dummy\"},\n{\"eqIndex\":1,\"section\":\"initial\",\"tag\":\"assign\",\"defines\":[\"y\"],\"uses\":[\"time\",\"omega\",\"A\"],\"equation\":[\"A * sin(omega * time)\"],\"source\":{\"info\":{\"file\":\"<interactive>\",\"lineStart\":5,\"lineEnd\":5,\"colStart\":3,\"colEnd\":36},\"within\":[\"Real\"]}},\n{\"eqIndex\":2,\"section\":\"regular\",\"tag\":\"assign\",\"defines\":[\"$cse1\"],\"uses\":[\"time\",\"omega\"],\"equation\":[\"sin(omega * time)\"],\"source\":{\"info\":{\"file\":\"\",\"lineStart\":0,\"lineEnd\":0,\"colStart\":0,\"colEnd\":0}}},\n{\"eqIndex\":3,\"section\":\"regular\",\"tag\":\"assign\",\"defines\":[\"y\"],\"uses\":[\"$cse1\",\"A\"],\"equation\":[\"A * $cse1\"],\"source\":{\"info\":{\"file\":\"<interactive>\",\"lineStart\":5,\"lineEnd\":5,\"colStart\":3,\"colEnd\":36},\"within\":[\"Real\"]}}\n],\n\"functions\":[\n]\n}"